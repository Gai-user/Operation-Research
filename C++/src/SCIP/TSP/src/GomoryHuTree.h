/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   GomoryHuTree.h
 * @brief  Gomory-Hu 割树生成器 —— 为无向图中的全局割生成不等式
 * @author Georg Skorobohatyj
 * @author Timo Berthold
 *
 * @details
 * 本文件实现了 Gomory-Hu 树（Gomory-Hu Tree，亦称割树 / Cut Tree）的数据结构与算法。
 * Gomory-Hu 树是一棵带权树，它紧凑地表示了无向图中所有节点对之间的最小 s-t 割。
 * 在 TSP（旅行商问题）的 SCIP 求解过程中，该割树用于生成全局割平面（global cuts），
 * 以加强线性规划松弛。
 *
 * ## 核心组件
 *   - **GraphNode / GraphEdge**：无向图的邻接表表示，采用半边模型（half-edge）。
 *     每条无向边由两条方向相反的半边组成，支持残余网络操作。
 *   - **Graph**：图的容器，管理节点和边的生命周期。节点和半边分别存储在连续的
 *     std::vector 中，保证内存局部性和 cache 友好性。
 *   - **MaxFlowContext**：push-relabel 最大流算法实现。采用最高标签优先策略
 *     （highest-label）和全局重标记（global relabel）优化。
 *   - **ghc_tree()**：主入口，从带容量的图构建 Gomory-Hu 割树，提取分离割。
 *
 * ## 半边模型（Half-Edge Model）
 * ~~~
 *   无向边 (u, v)，容量 c
 *   ┌──────────────────────────────────┐
 *   │  半边 A: u → v, capacity = c     │
 *   │  半边 B: v → u, capacity = c     │
 *   │  A.reverse = &B,  B.reverse = &A │
 *   │  A.target = &v,   B.target = &u  │
 *   └──────────────────────────────────┘
 * ~~~
 * 在最大流算法中，capacity 表示正向还可推送的流量，residual_capacity 表示
 * 反向可推送的流量（即"撤消"已推送流量的能力）。
 *
 * ## 算法参考
 *   - Gomory, R. E. and Hu, T. C. (1961).
 *     "Multi-terminal network flows."
 *     Journal of the Society for Industrial and Applied Mathematics, 9(4), 551–570.
 *   - Gusfield, D. (1990).
 *     "Very simple methods for all pairs network flow analysis."
 *     SIAM Journal on Computing, 19(1), 143–155.
 *   - Goldberg, A. V. and Tarjan, R. E. (1988).
 *     "A new approach to the maximum-flow problem."
 *     Journal of the ACM, 35(4), 921–940.
 */

#pragma once

#include <vector>
#include <memory>
#include <queue>

#include "objscip/objscip.h"

/**
 * @brief 数值比较的容差值（epsilon）
 *
 * 在最大流算法和割的构造中，用于判断容量是否为零、残余容量是否存在等。
 * 所有涉及浮点容量比较的地方都使用此常量作为容差，以避免舍入误差的影响。
 *
 * @note 该值为 1.0E-10，在 double 精度下足够小，同时能稳定处理 TSP 中常见的边长数值。
 */
constexpr double GH_EPS = 1.0E-10;


/** 前向声明，GraphNode 需要引用 GraphEdge */
struct GraphEdge;

/**
 * @brief 无向图的节点结构体
 *
 * 表示图中的一个顶点。该结构体服务于两个关键角色：
 *   1. **静态属性**：存储节点的 ID、坐标（用于 TSP 可视化）
 *   2. **算法运行时数据**：作为 push-relabel 最大流算法中的工作节点，
 *      维护距离标签（dist_label）、超额流（excess_flow）、活跃状态（is_active）等
 *
 * 节点通过 first_edge 指针维护其邻接表（边的单向链表）。
 * 在 Gomory-Hu 树构造过程中，parent 指针用于维护树形结构。
 */
struct GraphNode
{
   /// 节点编号（0-based，唯一标识）
   int  id                  = 0;

   /**
    * @brief 距离标签（distance label）
    *
    * push-relabel 算法中用于引导推流方向的核心数据结构。
    * 流量只能从距离标签较高的节点推向距离标签较低的节点（沿"下坡"方向）。
    * 汇点的距离标签始终为 0，源点的初始距离标签为节点总数 n。
    * 当节点无法推送时（所有出边饱和），通过 relabel 操作提升其距离标签。
    */
   int  dist_label          = 0;

   /// TSP 平面坐标 X（用于计算欧几里得距离和可视化）
   double coord_x           = 0.0;
   /// TSP 平面坐标 Y
   double coord_y           = 0.0;

   /**
    * @brief 超额流（excess flow）
    *
    * push-relabel 算法的核心概念。当流入节点的流量大于流出时，差额在此累积。
    * 算法维护"预流"（preflow）而非"流"（flow）——预流允许节点有超额流，
    * 而传统流要求每个中间节点的流入等于流出。算法终止时所有中间节点
    * 超额流归零，预流变为合法流。
    */
   double excess_flow       = 0.0;

   /**
    * @brief 最小割容量（minimum cut capacity）
    *
    * 在 Gomory-Hu 树的构建过程中，记录从该节点到其 parent 节点的最小 s-t 割容量。
    * 这是 Gomory-Hu 树边的权重，表示分离该节点与父节点的最小代价。
    */
   double min_cut_capacity  = 0.0;

   /**
    * @brief 可达性标记
    *
    * 在最大流算法的 BFS 阶段，从汇点出发沿残余网络反向搜索，
    * 可达的节点标记为 true。最终 is_reachable == true 的节点在汇侧，
    * is_reachable == false 的节点在源侧，两者之间形成最小割。
    * 初始值为 false（不可达，假设在源侧）。
    */
   bool  is_reachable       = false;

   /**
    * @brief 活跃标志
    *
    * push-relabel 算法中标记节点是否仍在参与计算。
    * 当节点因 gap heuristic 被判定为与汇点不可达时，is_active 设为 false，
    * 后续迭代将跳过该节点，这是重要的性能优化。
    */
   bool  is_active          = true;

   /// 邻接链表头指针，指向节点的第一条出边
   GraphEdge* first_edge    = nullptr;

   /**
    * @brief 当前扫描边指针
    *
    * discharge 操作中用于从上一次停止的位置继续扫描邻接边。
    * 避免每次 discharge 都从头遍历邻接表，是 push-relabel 算法的重要优化。
    */
   GraphEdge* current_edge  = nullptr;

   /**
    * @brief 活跃节点栈的 next 指针
    *
    * push-relabel 算法将活跃节点按距离标签分层组织在桶（bucket）中，
    * 同一层的节点通过 stack_next 连接成单向链表，实现 LIFO 栈语义。
    */
   GraphNode* stack_next    = nullptr;

   /**
    * @brief 父节点指针
    *
    * 在 Gomory-Hu 树构造过程中，记录该节点在树中的父节点。
    * 初始时所有节点（除根节点 0 外）的 parent 指向节点 0。
    * 随着算法进行，parent 指针会被更新以构建割树。
    */
   GraphNode* parent        = nullptr;

   /**
    * @brief 判断节点是否有出边（即是否连接到图中）
    * @return true 如果节点至少有一条邻接边，否则 false
    */
   bool hasEdge() const noexcept { return first_edge != nullptr; }

   /**
    * @brief 在邻接边中查找第一条满足谓词的边
    * @tparam Pred 谓词类型，接受 GraphEdge* 返回 bool
    * @param pred  谓词函数对象
    * @return 第一条满足 pred(e) == true 的边，若不存在则返回 nullptr
    *
    * 典型用法：查找容量大于某阈值的边、查找指向特定节点的边等。
    */
   template<typename Pred>
   GraphEdge* findEdge(Pred pred) const;

   /**
    * @brief 将 current_edge 向前推进，跳过容量 <= eps 的边
    * @param eps 容量阈值，默认为 GH_EPS
    *
    * 在 discharge 操作中，需要找到第一条"可推送"的边（即容量为正）。
    * 此函数跳过所有饱和边（capacity <= eps），加速后续的扫描过程。
    * 如果整条链表均饱和，current_edge 被置为 nullptr，表示需要执行 relabel。
    */
   void advanceCurrentEdge(double eps = GH_EPS) noexcept;

   /**
    * @brief 遍历所有容量大于阈值 eps 的出边
    * @tparam Func 回调函数类型，签名为 void(GraphEdge*)
    * @param func  对每条满足条件的边调用的回调函数
    * @param eps   容量阈值，默认 GH_EPS
    *
    * 对邻接链表中的每条出边，如果其 capacity > eps，则调用 func(edge)。
    * 可用于统计、收集、或批量操作满足容量条件的边。
    *
    * @note 不会修改邻接表结构，回调函数只读或通过指针修改边属性。
    *
    * 典型用法：
    * @code
    *   node.forEachOutgoingEdge([&](GraphEdge* e) {
    *       total_capacity += e->capacity;
    *   });
    * @endcode
    */
   template<typename Func>
   void forEachOutgoingEdge(Func&& func, double eps = GH_EPS) const
      noexcept(noexcept(func(std::declval<GraphEdge*>())));

   template<typename Func>
   void forEachIncomingEdge(Func&& func, double eps = GH_EPS) const
      noexcept(noexcept(func(std::declval<GraphEdge*>())));
};


/**
 * @brief 无向图的边结构体（半边 / half-edge）
 *
 * 每条无向边由两条"半边"表示，两条半边互为反向（通过 reverse 指针关联）。
 * 半边中存储了容量、残余容量和边长信息。
 *
 * ## 容量与残余容量
 * 以无向边 (u, v)，原始容量 C 为例：
 *   - **半边 A** (u→v): capacity = C - flow,  residual_capacity = flow
 *   - **半边 B** (v→u): capacity = C - flow,  residual_capacity = flow
 *
 * 初始时 flow = 0，两条半边的 capacity = C, residual_capacity = 0。
 * 当从 u 推送 f 单位流量到 v 时：
 *   - A.capacity 减少 f，A.residual_capacity 增加 f（允许"撤回"）
 *   - B.capacity 增加 f（实际是用 B.residual_capacity 来推），B.residual_capacity 减少 f
 *
 * ## 邻接链表
 * 每条半边通过 next 指针形成节点邻接表的单向链表，
 * 通过 reverse 指针指向其反向半边，通过 target 指针指向目标节点。
 */
struct GraphEdge
{
   /**
    * @brief 正向容量（capacity）
    *
    * 沿该边方向还可推送的最大流量。在残余网络中，capacity > 0 表示
    * 该方向上仍有可用容量。初始值为无向边的总容量。
    */
   double capacity          = 0.0;

   /**
    * @brief 残余容量（residual capacity）
    *
    * 沿该边反向可推送的最大流量（即已推送流量的大小）。
    * 在残余网络中，residual_capacity > 0 表示存在反向弧，
    * 允许流量"回流"或"撤回"。这是最大流算法能够纠正"错误"推送的关键机制。
    */
   double residual_capacity = 0.0;

   /// 边长（在 TSP 中为两节点之间的欧几里得距离或其他度量）
   double length            = 0.0;

   /// 同一节点邻接链表中的下一条边（单向链表）
   GraphEdge* next          = nullptr;

   /**
    * @brief 反向半边指针
    *
    * 每条无向边由一对互为反向的半边组成。reverse 指针指向配对的另一条半边。
    * 两条半边总是成对出现（如半边 e 和 e+num_edges）。
    */
   GraphEdge* reverse       = nullptr;

   /**
    * @brief 目标节点指针
    *
    * 指向该边的目标节点（即流量流向的终点）。
    * 在无向图中，另一条半边的 target 指向相反方向的节点。
    */
   GraphNode* target        = nullptr;

   /**
    * @brief 关联的 SCIP 决策变量
    *
    * 每条无向边对应一个 SCIP 二元变量 x_e ∈ {0, 1}，
    * 表示该边是否被选入 TSP 环。两条半边共享同一个变量引用。
    */
   SCIP_VAR* scip_var       = nullptr;

   /**
    * @brief 判断该边方向是否具有有效的正向容量（> eps）
    * @param eps 数值容差，默认 GH_EPS
    * @return true 如果 capacity > eps，即该方向上仍有可推送的流量
    */
   bool hasCapacity(double eps = GH_EPS) const noexcept { return capacity > eps; }

   /**
    * @brief 判断该边反向是否具有残余容量（> eps）
    * @param eps 数值容差，默认 GH_EPS
    * @return true 如果 residual_capacity > eps，即可以"撤消"该方向上的流量
    *
    * 在残余网络中，residual_capacity > 0 表示存在反向弧，允许流量回流。
    * 这是最大流算法能够纠正"错误"推送的关键机制。
    */
   bool hasResidualCapacity(double eps = GH_EPS) const noexcept { return residual_capacity > eps; }

   /**
    * @brief 设置两条互为反向的半边引用关系
    * @param other 反向半边指针
    *
    * 同时设置本边的 reverse 指向 other，以及 other 的 reverse 指向本边。
    * 调用后，两条半边形成一个完整的无向边对。
    *
    * @pre other != nullptr，且 other 尚未与其他半边配对
    */
   void setReverseEdge(GraphEdge* other) noexcept
   {
      reverse = other;
      other->reverse = this;
   }

   /**
    * @brief 将本边插入到节点的邻接链表头部
    * @param node 目标节点
    *
    * 采用头插法（O(1)），不维护任何顺序。邻接链表为单向链表，
    * 通过本函数构建的链表可通过 next 指针遍历。
    *
    * @pre node != nullptr
    */
   void linkTo(GraphNode* node) noexcept
   {
      next = node->first_edge;
      node->first_edge = this;
   }

   /**
    * @brief 将 SCIP 变量关联到本边及其反向边
    * @param scip SCIP 实例指针
    * @param v    要关联的 SCIP 变量
    * @return SCIP_OKAY 如果操作成功，否则返回 SCIP 错误码
    *
    * 一条无向边关联一个 SCIP 变量，但在数据模型中两个半边都持有引用。
    * 该函数同时对本边和反向边的 scip_var 指针调用 SCIPcaptureVar 来增加引用计数，
    * 防止变量在边仍使用时被 SCIP 释放。
    *
    * @pre scip != nullptr, v != nullptr, reverse != nullptr
    */
   SCIP_RETCODE setVar(SCIP* scip, SCIP_VAR* v)
   {
      assert(scip != nullptr);
      assert(v != nullptr);
      assert(reverse != nullptr);

      scip_var = v;
      SCIP_CALL(SCIPcaptureVar(scip, scip_var));

      reverse->scip_var = v;
      SCIP_CALL(SCIPcaptureVar(scip, reverse->scip_var));

      return SCIP_OKAY;
   }

   /**
    * @brief 释放本边及其反向边持有的 SCIP 变量引用
    * @param scip SCIP 实例指针
    * @return SCIP_OKAY 如果操作成功，否则返回 SCIP 错误码
    *
    * 减记变量的引用计数。先释放反向边再释放本边，遵循与 setVar 对称的释放顺序。
    * 在图的析构或边删除时调用，释放 SCIP 资源。
    *
    * @pre scip != nullptr, reverse != nullptr
    */
   SCIP_RETCODE releaseVar(SCIP* scip)
   {
      assert(scip != nullptr);
      assert(reverse != nullptr);

      SCIP_CALL(SCIPreleaseVar(scip, &reverse->scip_var));
      SCIP_CALL(SCIPreleaseVar(scip, &scip_var));

      return SCIP_OKAY;
   }
};


/* ---- GraphNode 内联成员（需在 GraphEdge 定义之后实现）---- */

/**
 * @brief 遍历邻接边，返回第一条满足谓词的边
 * @tparam Pred 谓词函数对象类型，签名为 bool(GraphEdge*)
 * @param pred  谓词实例
 * @return 满足 pred(e) == true 的第一条边，若所有边都不满足则返回 nullptr
 *
 * 线性遍历节点的邻接表。典型的谓词包括：
 *   - 查找容量大于 eps 的边
 *   - 查找指向特定目标节点的边
 *   - 查找未饱和的残余边
 */
template<typename Pred>
GraphEdge* GraphNode::findEdge(Pred pred) const
{
   for (GraphEdge* e = first_edge; e != nullptr; e = e->next)
      if (pred(e))
         return e;
   return nullptr;
}

/**
 * @brief 将 current_edge 推进到第一条容量 > eps 的边
 * @param eps 容量阈值，默认 GH_EPS
 *
 * 推进逻辑：
 * ~~~
 *   while (current_edge != nullptr && current_edge->capacity <= eps)
 *      current_edge = current_edge->next;
 * ~~~
 *
 * 在 push-relabel 算法中，discharge 操作需要找到可推送的邻接边。
 * 此函数预先跳过所有饱和边，使 current_edge 始终指向第一条候选推送边。
 * 如果整条链表均饱和，则 current_edge 被置为 nullptr，表示当前节点无可用出边，
 * 需要执行 relabel 操作。
 */
inline void GraphNode::advanceCurrentEdge(double eps) noexcept
{
   /* 使用 capacity 而非 residual_capacity，与 C 版本中的 cap 语义保持一致。
    * capacity 在 maxflow 过程中是常量（不会被 push 修改），
    * 而 residual_capacity 会随流量推送而变化。
    * 此函数的作用是跳过"始终为零容量"的边，而非跳过"当前已饱和"的边。 */
   while (current_edge != nullptr && current_edge->capacity <= eps)
      current_edge = current_edge->next;
}


/**
 * @brief 遍历所有容量大于阈值 eps 的出边
 * @tparam Func 回调函数类型，签名为 void(GraphEdge*)
 * @param func  对每条满足条件的边调用的回调函数
 * @param eps   容量阈值，默认 GH_EPS
 *
 * 对邻接链表中的每条出边，如果其 capacity > eps，则调用 func(edge)。
 * 可用于统计、收集、或批量操作满足容量条件的边。
 *
 * @note 不会修改邻接表结构，回调函数只读或通过指针修改边属性。
 */
template<typename Func>
void GraphNode::forEachOutgoingEdge(Func&& func, double eps) const
   noexcept(noexcept(func(std::declval<GraphEdge*>())))
{
   for (GraphEdge* e = first_edge; e != nullptr; e = e->next)
      if (e->capacity > eps)
         func(e);
}


/**
 * @brief 遍历所有具有残余容量的入边
 * @tparam Func 回调函数类型，签名为 void(GraphEdge*)
 * @param func  对每条满足条件的边调用的回调函数
 * @param eps   容量阈值，默认 GH_EPS
 *
 * 在残余网络中，"入边"的方向判断通过反向半边的残余容量来实现。
 * 对邻接链表中的每条边 e，如果 e->reverse->residual_capacity > eps，
 * 表示可以从 e->target 推送流量"撤回"到当前节点。
 *
 * 这是 push-relabel 算法 BFS 阶段的核心操作：从汇点沿残余网络反向搜索。
 *
 * @note 该函数遍历的是当前节点的出边链表，但检查的是反向边的残余容量。
 */
template<typename Func>
void GraphNode::forEachIncomingEdge(Func&& func, double eps) const
   noexcept(noexcept(func(std::declval<GraphEdge*>())))
{
   for (GraphEdge* e = first_edge; e != nullptr; e = e->next)
      if (e->reverse->residual_capacity > eps)
         func(e);
}


/**
 * @brief 无向图容器类
 *
 * 管理图的所有节点和边。每条无向边由两条半边（half-edge）组成，
 * 因此存储 m 条半边即表示 m/2 条无向边。
 *
 * ## 设计约束
 *   - **不可拷贝**（= delete）：防止指针引用关系被浅拷贝破坏
 *   - **可移动**（= default）：支持高效的传递语义
 *   - 通过工厂方法 Graph::create() 创建实例，禁止直接构造
 *
 * ## 内部存储
 *   - nodes_：节点的连续数组（std::vector<GraphNode>），索引 = 节点 ID，保证内存局部性
 *   - edges_：半边的连续数组（std::vector<GraphEdge>），半边成对分配：
 *     无向边 k 由半边 2k 和 2k+1 组成（或半边 k 和 k+num_edges，取决于实现约定）
 *
 * ## 内存布局示意
 * ~~~
 *   nodes_: [Node0][Node1][Node2]...[NodeN-1]
 *             ↓      ↓      ↓
 *   edges_: [0→1][0→2][1→0][1→2][2→0][2→1]...
 *            半边0  半边1  半边2  半边3  半边4  半边5
 *             ←无向边0→  ←无向边1→  ←无向边2→
 * ~~~
 */
class Graph
{
public:
   /**
    * @brief 工厂方法：分配并构造包含 n 个节点和 m 条半边（即 m/2 条无向边）的图
    * @param num_nodes      节点数量
    * @param num_half_edges 半边数量（必须为偶数，每条无向边贡献两条半边）
    * @return 图的 unique_ptr，由调用方持有所有权
    *
    * 使用 unique_ptr 返回确保所有权清晰，避免内存泄漏。
    * 构造函数为 private，只能通过此工厂方法创建实例。
    */
   static std::unique_ptr<Graph> create(int num_nodes, int num_half_edges)
   {
      return std::unique_ptr<Graph>(new Graph(num_nodes, num_half_edges));
   }

   /// 禁止拷贝 —— 图内部包含指针交叉引用，浅拷贝会导致悬挂指针
   Graph(const Graph&) = delete;
   Graph& operator=(const Graph&) = delete;
   /// 允许移动 —— 支持按值返回和移动语义
   Graph(Graph&&) = default;
   Graph& operator=(Graph&&) = default;

   // --- 只读访问器 ---

   /** @return 节点数量 */
   int numNodes()        const noexcept { return nnodes_; }
   /** @return 无向边数量（半边数的一半） */
   int numEdges()        const noexcept { return nedges_; }
   /** @return 非零容量边的数量 */
   int numEdgesNonZero() const noexcept { return nedgesnonzero_; }
   /** @return 节点数组的 const 指针（线性存储，可索引访问） */
   const GraphNode* nodes() const noexcept { return nodes_.data(); }
   /** @return 半边数组的 const 指针 */
   const GraphEdge* edges() const noexcept { return edges_.data(); }
   /** @return 第 i 个节点的 const 引用 */
   const GraphNode& node(int i) const noexcept { return nodes_[i]; }
   /** @return 第 i 条半边的 const 引用 */
   const GraphEdge& edge(int i) const noexcept { return edges_[i]; }

   // --- 非 const 访问（供算法内部修改图状态）---

   /** @return 节点数组的可变指针 */
   GraphNode* nodes() noexcept { return nodes_.data(); }
   /** @return 半边数组的可变指针 */
   GraphEdge* edges() noexcept { return edges_.data(); }
   /** @return 第 i 个节点的可变引用 */
   GraphNode& node(int i) noexcept { return nodes_[i]; }
   /** @return 第 i 条半边的可变引用 */
   GraphEdge& edge(int i) noexcept { return edges_[i]; }

   /** @brief 设置非零容量边的计数 */
   void setNumEdgesNonZero(int n) noexcept { nedgesnonzero_ = n; }
   /** @return 非零容量边计数的引用（允许算法内部修改） */
   int&  numEdgesNonZeroRef()  noexcept { return nedgesnonzero_; }
   /** @return 节点计数的引用 */
   int&  numNodesRef()         noexcept { return nnodes_; }
   /** @return 边计数的引用 */
   int&  numEdgesRef()         noexcept { return nedges_; }

private:
   /**
    * @brief 私有构造函数 —— 仅由工厂方法 create() 调用
    * @param n 节点数量，用于预分配 nodes_ vector
    * @param m 半边数量，用于预分配 edges_ vector
    *
    * 初始化后所有节点的字段为零值 / 默认值。调用方负责后续填充坐标、
    * 连接邻接表等工作。
    */
   Graph(int n, int m)
      : nodes_(n), edges_(m), nnodes_(n), nedges_(m / 2), nedgesnonzero_(m / 2)
   {}

   /// 节点数组，索引为节点 ID
   std::vector<GraphNode> nodes_;
   /// 半边数组，下标 i 和 i+num_edges 构成一对互为反向的半边
   std::vector<GraphEdge> edges_;
   /// 节点数量缓存（等于 nodes_.size()）
   int nnodes_;
   /// 无向边数量缓存（等于 edges_.size() / 2）
   int nedges_;
   /// 容量非零的无向边数量（等于 nedges_ 减去零容量边的数量）
   int nedgesnonzero_;
};


/**
 * @brief push-relabel 最大流算法上下文
 *
 * 实现了用于计算 s-t 最大流的 push-relabel 算法（也称为 preflow-push 算法）。
 * 与经典的 Ford-Fulkerson / Edmonds-Karp 算法相比，push-relabel 在实践中
 * 通常更快，尤其适用于稠密图。在最坏情况下时间复杂度为 O(V²√E)。
 *
 * ## 算法核心思想
 *
 * 传统增广路算法（Ford-Fulkerson）每次迭代寻找一条从源到汇的路径并推送流量，
 * 每次需要完整的图遍历。
 *
 * Push-relabel 采用"局部视角"：每个节点只关心其直接邻接节点，通过以下操作
 * 逐步将超额流推向汇点：
 *
 *   1. **预流阶段（preflow）**：将源点的所有出边容量一次性推送出去，
 *      使源点的邻接节点获得超额流。此时图的"流"不满足守恒律（中间节点
 *      流入 ≠ 流出），称为"预流"（preflow）。
 *
 *   2. **推进-重标记阶段（push-relabel）**：
 *      - **push**：活跃节点（excess_flow > 0）将超额流沿可行边推送至
 *        距离标签更低的邻接节点。推送量 = min(超额流, 边容量)。
 *      - **relabel**：当活跃节点的所有出边都饱和（无法推送），提升其
 *        距离标签为所有邻接节点中最小 dist_label + 1，以"解除阻塞"。
 *
 *   3. **算法终止**：所有节点（除源和汇外）的超额流为零时，
 *      预流变为合法流，汇点的超额流即为最大流值。此时从源点出发沿残余
 *      网络 BFS 可达的节点集合构成最小割的源侧。
 *
 * ## 关键优化
 *
 *   - **最高标签优先策略（highest-label）**：总是处理 dist_label 最大的活跃节点。
 *     将活跃节点按距离标签分层存储在桶数组中，每次 O(1) 取出最高层节点。
 *     避免了低效的 FIFO 或随机选择。
 *
 *   - **全局重标记（global relabel）**：周期性从汇点做反向 BFS，将所有可达节点
 *     的距离标签更新为到汇点的精确距离。这比逐步 relabel 收敛更快。
 *
 *   - **间隙启发式（gap heuristic）**：如果某层距离标签 k 的节点数为零（间隙），
 *     则所有 dist_label > k 且 < bound_ 的节点被标记为不可达（is_active = false），
 *     因为它们永远无法到达汇点。大幅减少无用操作。
 */
class MaxFlowContext
{
public:
   /**
    * @brief 构造最大流上下文
    * @param n 图中节点数
    *
    * active_ 和 number_ 分配为 n+1 大小以支持距离标签范围 [0, n]。
    * bound_ 初始为 n，co_check_ 初始为 true（首轮执行连通性检查）。
    */
   explicit MaxFlowContext(long n)
      : active_(n + 1, nullptr)
      , number_(n + 1, 0L)
      , max_dist_(0)
      , bound_(n)
      , co_check_(true)
   {}

   /**
    * @brief 计算从源点到汇点的最大流
    * @param graph  图指针（图结构中包含容量、坐标等信息）
    * @param source 源点指针
    * @param sink   汇点指针
    * @return 最大流量值（double）；若图不连通则返回 -1.0
    *
    * 主入口，完成以下步骤：
    *   1. 重置所有节点状态，恢复边的原始容量
    *   2. 连通性检查：若存在孤立节点，返回 -1.0
    *   3. 从汇点做反向 BFS，为每个节点建立精确距离标签
    *   4. 从源点推送预流（饱和所有源点出边）
    *   5. 迭代执行 discharge 操作直到活跃节点列表为空
    *   6. 返回汇点的超额流 = 最大流值
    *
    * 调用后，节点的 is_reachable 标记可用于构造最小割：
    * is_reachable == true 的节点在汇侧，is_reachable == false 的在源侧。
    */
   double maxflow(Graph* graph, GraphNode* source, GraphNode* sink);

   /**
    * @brief 执行全局重标记（global relabel）
    * @param graph 图指针
    * @param sink  汇点指针
    *
    * 从汇点做反向 BFS，将每个可达节点的 dist_label 重置为到汇点的精确距离。
    * 不可达节点的 dist_label 重置为 n（节点总数）。
    * 此操作周期性执行，用于收紧距离标签，避免因逐步 relabel 导致标签膨胀。
    */
   void globalRelabel(Graph* graph, GraphNode* sink);

private:
   /**
    * @brief 从汇点执行 BFS，建立初始精确距离标签
    * @param graph 图指针
    * @param sink  汇点指针
    *
    * 初始化阶段调用。以汇点为起点沿残余网络做反向 BFS（通过检查反向边的
    * residual_capacity > 0），为每个可达节点设置 dist_label = BFS 层数。
    * 不可达节点的 dist_label 保持为 n。
    *
    * BFS 过程中同时标记 is_reachable = true（可达节点在汇侧），
    * 并将已有超额流的可达节点推入活跃栈。
    */
   void bfsFromSink(Graph* graph, GraphNode* sink);

   /**
    * @brief 从源点推送初始预流
    * @param graph  图指针
    * @param source 源点指针
    * @param sink   汇点指针
    *
    * 初始化阶段调用。将源点到邻接节点的全部容量作为预流一次性推送出去：
    *   - 源点 dist_label 设为 n（保证不会被 push back）
    *   - 对源点的每条出边（capacity > 0），推送 capacity 单位的流量
    *   - 邻接节点获得对应的超额流，成为活跃节点
    *
    * 推送后源点被标记为 is_active = false，不再参与后续 discharge 循环。
    */
   void pushPreflow(Graph* graph, GraphNode* source, GraphNode* sink);

   /**
    * @brief 对活跃节点执行 discharge 操作
    * @param graph       图指针
    * @param active_node 活跃节点指针
    * @return true 如果节点仍活跃（未完成放电，需要 relabel），
    *         false 如果超额流已清零（放电完成）
    *
    * Discharge 是 push-relabel 的核心操作：
    * 反复尝试从 active_node 推送超额流到 dist_label = active_node.dist_label - 1
    * 的邻接节点。若 push 成功，检查邻接节点是否变为活跃并推入活跃栈。
    * 若所有出边饱和无法推送，返回 true 表示需要 relabel。
    *
    * 使用 current_edge 指针避免每次从头扫描邻接表。
    */
   bool discharge(Graph* graph, GraphNode* active_node);

   /**
    * @brief 沿边推送流量
    * @param from   源节点（推送方）
    * @param to     目标节点（接收方）
    * @param edge   推送沿用的边（from → to 方向）
    * @param amount 推送流量的大小
    *
    * 更新四个值：
    *   - from->excess_flow -= amount
    *   - to->excess_flow   += amount
    *   - edge->capacity    -= amount
    *   - edge->residual_capacity += amount  （允许撤回）
    *   - 反向边：capacity += amount, residual_capacity -= amount
    */
   void push(GraphNode* from, GraphNode* to, GraphEdge* edge, double amount);

   /**
    * @brief 对节点执行 relabel 操作（提升距离标签）
    * @param graph       图指针
    * @param active_node 待重标记的节点
    * @return true 如果重标记成功（节点仍活跃），
    *         false 如果节点因间隙启发式被标记为不活跃
    *
    * 将 active_node 的 dist_label 提升为其邻接节点（出边 capacity > 0）中
    * 最小 dist_label + 1。如果当前层在重标记后为空，触发间隙启发式：
    * 将该层以上且低于 bound_ 的所有节点标记为 is_active = false，
    * 因为它们在当前残余网络中不可达汇点。
    */
   bool relabelNode(Graph* graph, GraphNode* active_node);

   /**
    * @brief 将节点推入活跃栈（按距离标签分层）
    * @param node 节点指针
    *
    * 在 active_[node->dist_label] 链表的头部插入 node。
    * 同时更新 max_dist_（当前最大活跃距离层）。
    */
   void pushActive(GraphNode* node) noexcept;

   /**
    * @brief 从最大距离层弹出活跃节点
    * @return 弹出的节点指针，若所有层均为空则返回 nullptr
    *
    * 实现最高标签优先策略 —— 总是返回 dist_label 最大的活跃节点。
    * 如果当前最大层已空，递减 max_dist_ 直到找到非空层。
    */
   GraphNode* popActive() noexcept;

   /**
    * @brief 活跃节点分层桶
    *
    * active_[d] 指向距离标签为 d 的活跃节点链表头。
    * 各节点通过 stack_next 链接。使用数组实现分层桶（bucket），
    * 支持 O(1) 入栈和取出最高层节点。
    */
   std::vector<GraphNode*> active_;

   /**
    * @brief 各距离层的节点计数
    *
    * number_[d] 记录距离标签为 d 的节点数量（包括非活跃节点）。
    * 在 relabel 操作和 global relabel 中用于维护层级信息和检测间隙。
    */
   std::vector<long>       number_;

   /// 当前活跃节点栈中的最大距离层（用于 highest-label 选择）
   long max_dist_;

   /**
    * @brief 距离上界
    *
    * bound_ 的值为节点总数 n。dist_label >= bound_ 的节点被认为与汇点不可达。
    * 间隙启发式中使用此值：如果某层在 [0, bound_) 区间为空，
    * 则将所有 dist_label 高于该层且低于 bound_ 的节点标记为 killed。
    */
   long bound_;

   /**
    * @brief 连通性检查标志
    *
    * 当 co_check_ 为 true 时，下一次 maxflow() 调用将执行连通性检查。
    * 首轮检查后设为 false，避免重复检查。
    */
   bool co_check_;

   /**
    * @brief discharge 操作执行计数器
    *
    * 用于调度 global relabel：每隔一定次数的 discharge 后触发一次
    * global relabel，平衡重标记开销与收敛速度。
    */
   long discharge_count_ = 0;
};


/**
 * @brief Gomory-Hu 割树构建算法（Gusfield 简化版）
 *
 * 从给定的带容量无向图中计算 Gomory-Hu 割树，同时提取分离割作为 TSP 不等式。
 *
 * ## Gusfield 算法流程
 *
 * 朴素的 Gomory-Hu 算法需要计算 n(n-1)/2 次最大流（每对节点一次）。
 * Gusfield (1990) 发现只需要 n-1 次最大流即可构建完整的 Gomory-Hu 树：
 *
 *   1. **初始化**：所有节点（除节点 0 外）的 parent 指向节点 0
 *   2. **主循环**（i = 1, 2, ..., n-1）：
 *      a. 以节点 i 为源点，parent[i] 为汇点，计算 s-t 最大流
 *      b. 从源点沿残余网络 BFS，标记可达节点集合
 *      c. 对于 j > i 且 parent[j] == parent[i] 且 j 在源侧可达的节点，
 *         将 parent[j] 更新为 i（树的重连）
 *   3. 最终 parent 数组构成一棵树，树边的权重 = maxflow(i, parent[i])
 *
 * ## 为什么只需 n-1 次？
 *
 * Gomory-Hu 树的关键性质：对于任意三个节点 a, b, c，它们之间的三个割中，
 * 最小的两个相等。这意味着 n 个节点的所有割信息可以被一棵 n-1 条边的树
 * 完全编码。
 *
 * Gusfield 的观察是：用特定的节点对顺序计算最大流并更新 parent 指针，
 * 就能在不额外计算的情况下得到正确的树结构。
 *
 * @param graph        输入图（非 const，算法会修改容量以计算流，每次迭代后重置）
 * @param cuts         输出参数：割的位向量集合。
 *                     cuts[i][j] = true 表示第 i 个割将节点 j 分配到汇侧
 * @param num_cuts     输出参数：生成的割数量
 * @param min_violation 违反阈值 —— 仅当割的容量 < 2.0 - min_violation 时
 *                     才将其记录为有效割（即违反 LP 解的割）
 * @return true 如果至少有一个割满足违反阈值，否则返回 false
 *
 * @note 该函数通过修改 graph 中的边容量来计算流，并在每次迭代后重置。
 *       因此 graph 必须可修改，且原始容量数据需要在调用前通过其他方式保存。
 *
 * @see Gusfield, D. (1990). "Very simple methods for all pairs network flow analysis."
 *      SIAM Journal on Computing, 19(1), 143–155.
 */
bool ghc_tree(
   Graph*                        graph,
   std::vector<std::vector<bool>>& cuts,
   int*                          num_cuts,
   double                        min_violation
   );
