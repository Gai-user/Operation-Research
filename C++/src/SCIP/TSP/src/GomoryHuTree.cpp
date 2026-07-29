/**@file   GomoryHuTree.cpp
 * @brief  Gomory-Hu 割树生成器实现 —— push-relabel 最大流 + Gomory-Hu 树构造
 * @author Georg Skorobohatyj
 * @author Timo Berthold
 *
 * @details
 * 本文件包含两个核心算法的完整实现：
 *
 * ## 1. Push-Relabel 最大流算法 (MaxFlowContext)
 *
 * 算法维护"预流"（preflow）而非传统"流"。预流允许节点有超额流
 * （流入 > 流出），通过以下操作逐步将超额流推向汇点：
 *
 *   - **push**：活跃节点将超额流沿"下坡"边（dist_label 相差为 1）推送
 *   - **relabel**：当节点所有出边饱和时，提升其 dist_label 以解除阻塞
 *   - **全局重标记**：周期性从汇点 BFS，将 dist_label 更新为精确距离
 *   - **间隙启发式**：检测距离标签的"空隙"，标记空隙上方的节点为不可达
 *
 * 算法终止时，所有节点超额流归零，预流变为合法流，汇点的超额流 = 最大流值。
 *
 * ## 2. Gomory-Hu 割树 (ghc_tree)
 *
 * 用 Gusfield 简化算法计算 Gomory-Hu 树。从 n 个节点的完全图出发，
 * 执行 n-1 次 s-t 最大流计算，通过巧妙的 parent 指针重连，
 * 构建一棵能表示所有 O(n²) 对节点间最小割的树。
 *
 * 在 TSP 割平面生成中，如果某个割的容量小于 2.0（即 LP 解中跨越该割
 * 的变量值之和小于 2），则该割违反度约束，可以作为分离割添加到 SCIP 中。
 */

#include <cstdio>
#include <cassert>
#include <algorithm>
#include <memory>

#include "objscip/objscip.h"
#include "GomoryHuTree.h"

/* ==========================================================================
 * MaxFlowContext 实现
 *
 * 本节实现 push-relabel 最大流算法的各个组件。
 * 算法参考：Goldberg & Tarjan (1988), "A new approach to the maximum-flow
 * problem", Journal of the ACM, 35(4), 921-940.
 * ========================================================================== */

/**
 * @brief 将节点推入活跃栈（按距离标签分层）
 *
 * 使用头插法将节点插入对应 dist_label 层的链表头部。
 * 同时更新 max_dist_ 以维护最高标签优先策略。
 *
 * @param node 要推入活跃栈的节点
 */
void MaxFlowContext::pushActive(GraphNode *node) noexcept
{
    // 更新最大活跃距离层（highest-label 策略需要知道当前最大的 dist_label）
    if (node->dist_label > max_dist_)
        max_dist_ = node->dist_label;

    node->stack_next = active_[node->dist_label];
    active_[node->dist_label] = node;
}

/**
 * @brief 从汇点执行 BFS，建立初始精确距离标签
 *
 * 以汇点为起点，沿残余网络做反向 BFS。
 * 反向搜索的关键：检查每条出边的反向半边是否有残余容量
 * （即 e->reverse->residual_capacity > eps），这表示在残余网络中存在
 * 从 e->target 指向当前节点的弧。
 *
 * 边遍历边标记：
 *   - 将可达节点标记为 is_reachable = true（在汇侧）
 *   - 设置 dist_label = 到汇点的 BFS 层数
 *   - 如果节点已有超额流，推入活跃栈
 *
 * BFS 完成后，未标记的节点（is_reachable == false）在源侧。
 * 源侧与汇侧之间的边集合构成最小 s-t 割。
 *
 * @param graph 图指针
 * @param sink  汇点指针
 */
void MaxFlowContext::bfsFromSink(Graph *graph, GraphNode *sink)
{
    int num_nodes = graph->numNodes();
    std::queue<GraphNode *> q;

    // 汇点自身总是在汇侧（可达），无论调用方是否预置了此标记
    sink->is_reachable = true;
    q.push(sink);

    int count = 1;

    while (!q.empty())
    {
        GraphNode *current = q.front();
        q.pop();

        int level = current->dist_label + 1;

        // 遍历当前节点的所有邻接边，检查反向边是否有残余容量
        for (GraphEdge *edge = current->first_edge; edge != nullptr; edge = edge->next)
        {
            GraphNode *neighbor = edge->target;
            // 三个条件同时满足才可沿残余网络反向搜索：
            // 1. neighbor->is_active：邻接节点仍在参与计算
            // 2. neighbor->is_reachable == false：尚未被访问
            // 3. edge->reverse->hasResidualCapacity(GH_EPS)：残余网络存在反向弧
            if (neighbor->is_active && !neighbor->is_reachable && edge->reverse->hasResidualCapacity(GH_EPS))
            {
                // 标记为汇点可达
                neighbor->is_reachable = true;
                neighbor->dist_label = level;
                ++count;
                ++number_[level];

                // 若是该节点有超额流量，推入活跃栈
                if (neighbor->excess_flow > GH_EPS)
                {
                    pushActive(neighbor);
                    if (level > max_dist_)
                    {
                        max_dist_ = level;
                    }
                }
                q.push(neighbor);
            }
        }
    }
    // BFS 不可达的活跃节点：dist_label 设为 n，标记为不活跃
    // 这些节点在残余网络中已被隔断，无法到达汇点
    for (int i = 0; i < num_nodes; ++i)
    {
        auto &node = graph->node(i);
        if (!node.is_reachable && node.is_active)
        {
            node.dist_label = num_nodes;
            node.is_active = false;
        }
    }

    discharge_count_ = count;
    if (count < bound_)
    {
        bound_ = count;
    }
}

/**
 * @brief 全局重标记（global relabel）
 *
 * 从汇点沿残余网络做反向 BFS，将所有可达节点的 dist_label 更新为
 * 到汇点的精确距离。不可达节点（is_reachable == false 且 is_active == true）
 * 的 dist_label 被设为 n（节点总数），is_active 设为 false。
 *
 * 此操作周期性执行，将因逐步 relabel 而膨胀的距离标签"拉回"精确值，
 * 显著加速收敛。
 *
 * @param graph 图指针
 * @param sink  汇点指针
 */
void MaxFlowContext::globalRelabel(Graph *graph, GraphNode *sink)
{
    int num_nodes = graph->numNodes();

    // 重置所有节点， 标记为不可达（在源侧），清空活跃栈链接
    // 重置 current_edge 到邻接表头部
    for (int i = num_nodes - 1; i >= 0; --i)
    {
        auto &node = graph->node(i);
        node.is_reachable = false;
        node.stack_next = nullptr;
        node.current_edge = node.first_edge;
        node.advanceCurrentEdge(GH_EPS);
    }

    // 汇点可达自身
    sink->is_reachable = true;

    // 清空活跃栈的所有层和层节点计数
    for (int i = 0; i <= num_nodes; ++i)
    {
        active_[i] = nullptr;
        number_[i] = 0;
    }

    // 从汇点 BFS 建立精确距离标签
    bfsFromSink(graph, sink);

    // 不可达且活跃的节点：dist_label 设为 n，标记为不活跃
    for (int i = 0; i < num_nodes; ++i)
    {
        auto &node = graph->node(i);
        if (!node.is_reachable && node.is_active)
        {
            node.dist_label = num_nodes;
            node.is_active = false;
        }
    }
}

/**
 * @brief 沿边推送流量
 *
 * 从 from 节点推送 amount 单位的流量到 to 节点。
 * 更新四个值以维持残余网络的一致性：
 *   - from->excess_flow 减少（流出）
 *   - to->excess_flow   增加（流入）
 *   - edge->capacity    减少（正向容量被消耗）
 *   - edge->residual_capacity 增加（允许撤回，形成反向弧）
 *   - 反向边做相反的更新
 *
 * @param from   源节点（推送方）
 * @param to     目标节点（接收方）
 * @param edge   推送沿用的边（from → to 方向）
 * @param amount 推送流量的大小
 */
void MaxFlowContext::push(GraphNode *from, GraphNode *to, GraphEdge *edge, double amount)
{
    // 从推送方扣除超额流
    from->excess_flow -= amount;
    // 接收方获得超额流
    to->excess_flow += amount;

    // 正向边：容量减少，残余容量增加（允许撤回）
    edge->residual_capacity -= amount;
    // 反向边：容量增加（可以被利用），残余容量减少
    edge->reverse->residual_capacity += amount;
}

/**
 * @brief 从源点推送初始预流
 *
 * 将源点到所有邻接节点的全部容量作为预流一次性推送出去。
 * 这是 push-relabel 算法的初始化步骤：
 *   1. 将源点的 dist_label 设为 n（最高值，防止流量被推回源点）
 *   2. 对每条源点出边（capacity > 0），推送全部容量
 *   3. 邻接节点获得超额流，成为活跃节点
 *   4. 源点标记为 is_active = false（不再参与后续 discharge）
 *
 * 注意：推送后源点的 excess_flow 为负数（所有出边容量之和的相反数），
 * 但源点标记为不活跃后，这个负值不影响算法正确性。
 *
 * @param graph  图指针
 * @param source 源点指针
 * @param sink   汇点指针
 */
void MaxFlowContext::pushPreflow(Graph *graph, GraphNode *source, GraphNode *sink)
{
    max_dist_ = 0;

    // 遍历源点所有相邻的边
    for (GraphEdge *edge = source->first_edge; edge != nullptr; edge = edge->next)
    {
        // 只推送有容量的边
        if (edge->hasCapacity(GH_EPS))
        {
            double amount = edge->residual_capacity;
            double old_excess = edge->target->excess_flow;

            push(source, edge->target, edge, amount);
            // 如果邻接节点获得超额流且之前超额流 ≈ 0，推入活跃栈
            // 排除汇点（等价于 C 版本中的 nptr != t_ptr && nptr->excess <= cap + EPS）
            if (edge->target != sink && old_excess <= GH_EPS)
            {
                pushActive(edge->target);
            }
        }
    }
}

/**
 * @brief 从最大距离层弹出活跃节点（最高标签优先策略）
 *
 * 总是返回当前 dist_label 最大的活跃节点。
 * 如果当前最大层已空，递减 max_dist_ 直到找到非空层或归零。
 *
 * @return 弹出的节点指针，若所有层均为空则返回 nullptr
 */
GraphNode *MaxFlowContext::popActive() noexcept
{
    // 降级查找：当前最大层为空，递减直到找到非空层或归零
    while (max_dist_ > 0 && active_[max_dist_] == nullptr)
    {
        --max_dist_;
    }
    GraphNode *node = active_[max_dist_];
    if (node != nullptr)
    {
        active_[max_dist_] = node->stack_next;
    }
    return node;
}

/**
 * @brief 对活跃节点执行 discharge 操作
 *
 * Discharge 是 push-relabel 算法的核心循环。它反复尝试将活跃节点
 * 的超额流推送到邻接节点（沿"下坡"方向：dist_label 相差为 1）。
 *
 * 操作流程：
 *   1. 从 current_edge 开始扫描邻接表（避免每次从头扫描）
 *   2. 对每条 capacity > 0 的边，检查是否满足推送条件：
 *      active_node.dist_label == target.dist_label + 1
 *   3. 如果满足，推送 min(excess_flow, capacity) 单位的流量
 *   4. 推送后检查目标节点是否变为活跃
 *   5. 如果超额流清零 → 返回 false（放电完成）
 *   6. 如果扫描完所有边仍无法清零 → 返回 true（需要 relabel）
 *
 * 推送条件的物理意义：
 *   dist_label 代表节点到汇点的"估计距离"。流量只能沿估计距离
 *   减少 1 的方向推送（类似水往低处流）。这保证了算法朝正确方向
 *  （汇点方向）推进流量。
 *
 * @param graph       图指针
 * @param active_node 活跃节点指针
 * @return true 如果节点仍需要 relabel，false 如果放电完成
 */
bool MaxFlowContext::discharge(Graph *graph, GraphNode *active_node)
{
    // 从上次扫描的位置继续
    GraphEdge *edge = active_node->current_edge;
    while (edge != nullptr && active_node->excess_flow > GH_EPS)
    {
        GraphNode *target = edge->target;

        // 推送条件：
        // 1. 边有正向残余容量（edge->residual_capacity > 0）
        // 2. 距离标签差值为 1 （水向低处流）
        if (edge->residual_capacity > GH_EPS &&
            active_node->dist_label == target->dist_label + 1)
        {
            // 推送的流量 = min(节点超额流量，边的残余容量)
            double amount = std::min(active_node->excess_flow, edge->residual_capacity);
            double old_excess = target->excess_flow;
            push(active_node, target, edge, amount);

            // 目标节点获得流量后可能变为活跃节点
            // 只有在目标节点流量约为0时，将其推入活跃节点
            if (target->dist_label > 0 && old_excess <= GH_EPS)
            {
                pushActive(target);
            }
        }

        edge = edge->next;
    }

    // 更新扫描指针，记录下次从哪里开始
    active_node->current_edge = edge;

    // 超额流量仍未清零 -> 需要 relabel
    if (active_node->excess_flow > GH_EPS)
    {
        return true;
    }

    return false;
}

/**
 * @brief 对活跃节点执行 relabel 操作（提升距离标签）
 *
 * 当活跃节点的所有出边都饱和，无法继续推送时，需要"抬高"节点的
 * 距离标签以解除阻塞。
 *
 * Relabel 操作：
 *   1. 扫描所有出边（capacity > 0），找到邻接节点中最小 dist_label
 *   2. 将当前节点的 dist_label 设为 min_dist_label + 1
 *   3. 如果当前层变为空（间隙），触发间隙启发式
 *
 * 间隙启发式（Gap Heuristic）：
 *   如果某层 k 的节点数为零，意味着距离标签在 k 处出现了"空隙"。
 *   所有 dist_label > k 且 < bound_ 的节点的超额流永远无法到达汇点
 *   （因为流量只能沿 dist_label 逐层下降），将它们标记为 is_active = false。
 *   这是 push-relabel 最重要的优化之一。
 *
 * @param graph       图指针
 * @param active_node 待重标记的节点
 * @return true 如果重标记成功（节点仍活跃），false 如果节点被间隙启发式杀死
 */
bool MaxFlowContext::relabelNode(Graph *graph, GraphNode *active_node)
{
    int num_nodes = graph->numNodes();

    int min_label = num_nodes;
    int old_label = active_node->dist_label;

    active_node->current_edge = nullptr;
    for (GraphEdge *edge = active_node->first_edge; edge != nullptr; edge = edge->next)
    {
        if (edge->residual_capacity > GH_EPS)
        {
            int label = edge->target->dist_label;
            if (label < min_label)
            {
                min_label = label;

                if (active_node->current_edge == nullptr)
                {
                    active_node->current_edge = edge;
                }
            }
        }
    }

    int new_label = min_label + 1;

    // 边界检查：若 new_label >= bound_，则该节点与汇点不可达
    if (new_label >= bound_)
    {
        --number_[old_label];
        active_node->dist_label = num_nodes;
        active_node->is_active = false;
        --bound_;
        return false;
    }

    active_node->dist_label = new_label;

    --number_[old_label];
    ++number_[new_label];

    // 更新 max_dist_ 为新标签值
    // relabel 后当前节点的标签可能成为最高活跃标签
    if (new_label > max_dist_) 
    {
        max_dist_ = new_label;
    }

    // 不将节点推入活跃栈，而是继续处理这一节点
    // 间隙启发式：如果旧层变为空，标记更高层节点为不可达
    if (number_[old_label] == 0)
    {
        for (int i = 0; i < num_nodes; ++i)
        {
            auto& node = graph->node(i);
            int label = node.dist_label;

            if (label > old_label && label < bound_ && node.is_active)
            {
                node.is_active = false;
                --number_[label];;
            }
        }
        return false;
    }
    return true;
}

/**
 * @brief 计算从源点到汇点的最大流
 *
 * Push-relabel 算法的主入口。完整的最大流计算流程。
 *
 * @param graph  图（含节点和边容量，算法会修改残余容量）
 * @param source 源点
 * @param sink   汇点
 * @return 最大流量值（即汇点的超额流）；若图不连通则返回 -1.0
 *
 * ## 算法流程（逐步）
 *
 * ### 步骤 1：初始化
 * 重置所有节点的运行时状态：
 *   - current_edge 指向邻接表头部
 *   - excess_flow 清零
 *   - is_active = true（所有节点初始参与计算）
 *   - is_reachable = false（初始假设在源侧）
 *
 * ### 步骤 2：连通性检查
 * 如果存在孤立节点（first_edge == nullptr），图不连通，返回 -1.0。
 *
 * ### 步骤 3：重置边容量
 * 将所有半边的 residual_capacity 恢复为 capacity（即原始容量），
 * 确保每次 maxflow 调用从相同初始状态开始。
 *
 * ### 步骤 4：汇点 BFS
 * 从汇点沿残余网络做反向 BFS，建立精确距离标签。
 * 同时进行连通性检查：如果存在不可达节点，返回 -1.0。
 *
 * ### 步骤 5：源点预流推送
 * 将源点所有出边容量作为预流推送出去。源点标记为不活跃。
 *
 * ### 步骤 6：主循环
 * 反复执行以下操作直到活跃节点列表为空：
 *   a. 从最大距离层弹出一个活跃节点
 *   b. 对该节点执行 discharge 操作
 *   c. 如果 discharge 返回 true（需要 relabel），执行 relabel
 *   d. 周期性执行全局重标记
 *
 * ### 步骤 7：返回结果
 * 汇点的 excess_flow = 最大流值。
 * 节点的 is_reachable 标记可用于构造最小割。
 */
double MaxFlowContext::maxflow(Graph *graph, GraphNode *source, GraphNode *sink)
{
    int num_nodes = graph->numNodes();

    // ---- 步骤 1：重置边容量（必须在节点初始化之前，确保 advanceCurrentEdge 使用正确的 residual_capacity）----
    int num_half_edges = graph->numEdgesNonZero();
    int num_edges = graph->numEdges();

    // 将前半段半边（编号 0 ~ num_half_edges-1）的残余容量恢复为容量
    for (int i = 0; i < num_half_edges; ++i)
    {
        auto &edge = graph->edge(i);
        edge.residual_capacity = edge.capacity;
    }

    // 将后半段半边（编号 num_edges ~ num_edges+num_half_edges-1）同样处理
    for (int i = num_edges + num_half_edges - 1; i >= num_edges; --i)
    {
        auto &edge = graph->edge(i);
        edge.residual_capacity = edge.capacity;
    }

    // ---- 步骤 2：初始化所有节点 ----
    for (int i = 0; i < num_nodes; ++i)
    {
        auto &node = graph->node(i);

        // 重置扫描指针到邻接表头部，跳过始终为零容量的边
        node.current_edge = node.first_edge;
        node.advanceCurrentEdge(GH_EPS);

        // 孤立节点检测：没有任何邻接边的节点 → 图不连通
        if (node.current_edge == nullptr)
        {
            fprintf(stderr, "isolated node in input graph\n");
            return -1.0;
        }

        // 重置节点状态
        node.excess_flow = 0.0;
        node.stack_next = nullptr;
        node.is_active = true;
        node.is_reachable = false;
    }

    // ---- 步骤 3：清空活跃栈 ----
    for (int i = num_nodes; i >= 0; --i)
    {
        number_[i] = 0L;
        active_[i] = nullptr;
    }

    // ---- 步骤 4：汇点 BFS，建立初始距离标签 ----
    sink->dist_label = 0;
    sink->is_reachable = true;
    bfsFromSink(graph, sink);

    // 连通性检查：如果存在不可达节点（is_reachable == false 且 is_active == true），
    // 说明图不连通
    if (co_check_)
    {
        co_check_ = false;
        for (int i = 0; i < num_nodes; ++i)
        {
            auto &node = graph->node(i);
            if (!node.is_reachable && node.is_active)
            {
                return -1.0;
            }
        }
    }

    // ---- 步骤 5：源点预流推送 ----
    source->dist_label = num_nodes;
    sink->dist_label = 0;
    sink->excess_flow = 1.0;

    pushPreflow(graph, source, sink);

    source->is_active = false; // 源点不参与后续流程
    bound_ = num_nodes;
    discharge_count_ = 0;

        // 主循环
        while (max_dist_ > 0)
    {
        // 弹出标签最大的节点
        GraphNode *active_node = popActive();

        while (active_node != nullptr)
        {
            // 对当前活跃节点执行 discharge（扫描邻接边，推送超额流）
            // discharge 返回 true 表示超额流未清零，需要执行 relabel
            bool needs_relabel = discharge(graph, active_node);

            ++discharge_count_;

            // 周期性全局重标记：每 n 次 discharge 执行一次
            if (discharge_count_ >= num_nodes)
            {
                discharge_count_ = 0;
                globalRelabel(graph, sink);
            }

            // 如果 discharge 后仍有超额流量且节点活跃，执行relabel
            if (needs_relabel && active_node->is_active && active_node->excess_flow > GH_EPS)
            {
                // relabelNode 更新距离标签并重置 current_edge
                // 返回 true  → relabel 成功，节点仍活跃，current_edge 已重置，
                //               继续处理同一节点
                // 返回 false → 节点被间隙启发式杀死，处理下一节点
                if (relabelNode(graph, active_node))
                {
                    continue;
                }
            }

            active_node = popActive();
        }
    }

    return sink->excess_flow - 1.0;
}

/* ==========================================================================
 * Gomory-Hu 树构造
 *
 * 从带容量的无向图中构建 Gomory-Hu 割树并提取分离割。
 * 使用 Gusfield (1990) 的简化算法，只需 n-1 次最大流计算。
 * ========================================================================== */

/**
 * @brief 检查节点 v 是否在节点 root 到根节点的路径上
 *
 * 沿 parent 指针链向上追溯，判断节点 v 是否是 root 的祖先节点。
 * 在 Gomory-Hu 树中，这用于确定节点在树中的相对位置。
 *
 * @param graph 图指针
 * @param root  起点节点 ID
 * @param v     待检查节点 ID
 * @return true 如果 v 在 root 到根的路径上
 */
static bool nodeOnRootPath(Graph *graph, int i, int j)
{
    auto *nodes = graph->nodes();
    for (auto *node = &nodes[j]; node != &nodes[0]; node = node->parent)
    {
        if (node == &nodes[i])
        {
            return true;
        }
    }
    return false;
}

/**
 * @brief 从已构建的 Gomory-Hu 树中提取违反的割
 *
 * 遍历 GH 树的所有边（除了根节点到其父节点的虚拟边），
 * 对于每条边的 min_cut_capacity（即边的权重，等于最大流值），
 * 如果小于 2.0 - min_violation，说明该割违反 TSP 的度约束，
 * 将其记录为有效割。
 *
 * 割的构造方式：从节点 i 沿 parent 链向上找到根节点，路径上
 * 所有节点的 parent 方向定义了割的划分。
 *
 * @param graph        已构建 GH 树的图
 * @param cuts         输出参数：割的位向量集合
 * @param num_cuts     输出参数：有效割的数量
 * @param min_violation 违反阈值
 */
static void collectCuts(
    Graph *graph,
    std::vector<std::vector<bool>> &cuts,
    int *num_cuts,
    double min_violation)
{
    int k = 0;
    int num_nodes = graph->numNodes();

    cuts.resize(num_nodes - 1);

    // 检查每个非根节点的 min_cut_capacity
    for (int i = 1; i < num_nodes; ++i)
    {
        auto &node = graph->node(i);
        // 如果割容量 < 2.0 - min_violation，说明违反约束
        // 在 TSP 中，每个节点的度约束 = 2，所以割容量小于 2 意味着
        // LP 解中跨越该割的边变量值之和小于 2，是一个有效的分离割
        if (node.min_cut_capacity < 2.0 - min_violation)
        {
            cuts[k].resize(num_nodes);

            // 节点 0（根节点）始终在源侧（false）
            cuts[k][0] = false;

            // 对于其他节点，通过检查是否在节点 i 到根的路径上来划分
            for (int j = 1; j < num_nodes; ++j)
            {
                cuts[k][j] = nodeOnRootPath(graph, i, j);
            }
            ++k;
        }
    }
    cuts.resize(k);
    *num_cuts = k;
}

/**
 * @brief 构建 Gomory-Hu 割树并提取 TSP 分离割
 *
 * 使用 Gusfield (1990) 简化算法计算 Gomory-Hu 树。
 * 从 n 个节点的图出发，执行 n-1 次最大流计算。
 *
 * ## Gusfield 算法详解
 *
 * ### 为什么只需 n-1 次最大流？
 *
 * 考虑三个节点 a, b, c。定义 c(a,b) 为 a 和 b 之间的最小割容量。
 * Gomory-Hu 树的核心引理：c(a,b) ≥ min(c(a,c), c(b,c))。
 * 这意味着在三个节点中，三个割容量中最小的两个相等。
 * 因此，n 个节点之间的 O(n²) 个割信息可以被 n-1 条边的树完全编码。
 *
 * ### 算法步骤
 *
 * 1. **初始化**：所有非根节点（i=1..n-1）的 parent 指向节点 0（根）
 * 2. **主循环**（对 i = 1..n-1）：
 *    a. 以节点 i 为源点，parent[i] 为汇点，计算最大流
 *    b. 割容量 = 最大流值，存入 i.min_cut_capacity
 *    c. **树的重连**（核心步骤）：对于 j > i 且 parent[j] == 原 parent[i]
 *       且 j 在源侧（is_reachable == false）的节点，将 parent[j] 更新为 i
 *    d. 重置边容量，准备下一轮
 * 3. 从树中提取所有违反 TSP 约束的割
 *
 * ### 树的重连（Step 2c）的直观解释
 *
 * 最大流计算将图分为两部分：源侧（S）和汇侧（T）。
 * 对于 parent 指向汇点且在源侧的节点 j，节点 j 与节点 i 之间的割
 * 不大于 j 与汇点之间的割。因此，在 GH 树中，j 应该成为 i 的子节点
 * 而非原汇点的子节点。
 *
 * ### 异常处理
 *
 * 如果 maxflow 返回负值（图不连通），使用 BFS 标记构造一个降级割。
 * 虽然这会丢失 Gomory-Hu 树的结构信息，但至少能返回一个有效割。
 *
 * @param graph         输入图（会被修改）
 * @param cuts          输出：割的位向量集合
 * @param num_cuts      输出：有效割的数量
 * @param min_violation 违反阈值
 * @return true 如果至少有一个有效割
 */

bool ghc_tree(
    Graph *graph,
    std::vector<std::vector<bool>> &cuts,
    int *num_cuts,
    double min_violation)
{
    int num_nodes = graph->numNodes();

    // 创建最大流求解器
    MaxFlowContext maxflow_ctx(num_nodes);

    // 初始化：所有节点的parent指向根节点
    auto &root = graph->node(0);
    for (int i = num_nodes - 1; i >= 0; --i)
    {
        graph->node(i).parent = &root;
    }

    // gusfield 主循环：n - 1次最大流
    for (int i = 1; i < num_nodes; ++i)
    {
        GraphNode *source = &graph->node(i);
        GraphNode *sink = source->parent;

        // 计算source与sink之间的最大流
        double maxflow_value = maxflow_ctx.maxflow(graph, source, sink);

        // maxflow_value为负值表示图中存在孤立节点或不连通分量
        if (maxflow_value < 0.0)
        {
            // 降级策略：使用 BFS 留下的 is_reachable 标记构造一个割
            // is_reachable == true  → 汇侧（可从汇点沿残余网络到达）
            // is_reachable == false → 源侧（不可达）
            cuts.resize(1);
            cuts[0].resize(num_nodes);
            cuts[0][0] = false; // 根节点在源侧

            for (int j = 1; j < num_nodes; ++i)
            {
                GraphNode &node = graph->node(j);
                cuts[0][j] = !node.is_reachable;
            }
            *num_cuts = 1;
            return true;
        }

        // 更新该节点到父节点的最小割
        source->min_cut_capacity = maxflow_value;

        // ---- 树的重连（Gusfield 算法的核心）----
        // 对于所有 j ≠ i 的节点：
        //   如果 j 的 parent 等于 sink，且 j 在最大流后位于源侧
        //   （is_active == false，即不可从汇点到达），
        //   则将 j 的 parent 更新为 source（节点 i）
        //
        // 这维护了 GH 树的不变量：对于任意节点 j，
        // 割 (j, parent[j]) 的容量等于 j 与 parent[j] 之间的最小割。
        for (int j = 1; j < num_nodes; ++j)
        {
            auto &node = graph->node(j);
            if (&node != source && !node.is_active && node.parent == sink)
            {
                node.parent = source;
            }
        }

        // ---- GH 树旋转（Gomory-Hu 树的构造步骤）----
        // 如果 sink 的父节点在源侧（is_active == false），
        // 表示 sink 的父节点不可从汇点到达（与 sink 不在同一侧）。
        // 此时需要将 source 和 sink 在树中的位置进行旋转：
        //   1. source 成为 sink 父节点的子节点
        //   2. sink 成为 source 的子节点
        //   3. 交换两者的 min_cut_capacity
        //
        // 这确保 GH 树的不变量：
        //   树上每条边 (u, parent[u]) 的权重 = u 与 parent[u] 之间的最小割容量
        if (!sink->parent->is_active)
        {
            source->parent = sink->parent;
            sink->parent = source;
            source->min_cut_capacity = sink->min_cut_capacity;
            sink->min_cut_capacity = maxflow_value;
        }
    }

    // 收集所有违反的割
    collectCuts(graph, cuts, num_cuts, min_violation);

    return *num_cuts > 0;
}