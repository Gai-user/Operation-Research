/**@file   ProbDataTSP.cpp
 * @brief  TSP 问题数据类实现 —— 图的深拷贝与 SCIP 变量管理
 * @author Timo Berthold
 *
 * @details
 * 主要实现三个核心操作：
 *
 * ## 1. 图的深拷贝 (copyGraph)
 *
 * 创建 Graph 的完整副本，包括所有节点和半边。
 * 对于完全图，有 n*(n-1)/2 条无向边，即 n*(n-1) 条半边。
 * 半边在数组中成对存储：半边 k 和 k+num_edges 互为反向。
 *
 * ## 2. SCIP 变量管理
 *
 * 每条无向边关联一个 SCIP 二元变量。在图的拷贝和变换过程中，
 * 需要正确地将原始变量映射到新的/变换后的变量。
 * 使用 SCIPcaptureVar / SCIPreleaseVar 管理引用计数。
 *
 * ## 3. 资源释放
 *
 * 释放时先释放所有边关联的 SCIP 变量，再释放图本身。
 */

#include <cassert>
#include "objscip/objscip.h"
#include "ProbDataTSP.h"
#include "GomoryHuTree.h"

using namespace scip;

namespace tsp
{
    /**
     * @brief 创建图的深拷贝
     *
     * 将源图 src 的所有节点和边完整复制到新图 dest 中。
     * 复制内容：
     *   - 节点：ID、坐标、邻接表头（不复制 first_edge 指针，后续重建）
     *   - 半边：正向和反向半边的配对关系（setReverseEdge）、
     *          边长（length）、邻接表链接（linkTo）
     *
     * 半边存储布局（完全图）：
     * ~~~
     *   对每对节点 (i, j)，i < j：
     *     半边 k（正向 i→j）：存储在 edges[k]
     *     半边 k+num_edges（反向 j→i）：存储在 edges[k+num_edges]
     * ~~~
     *
     * @param dest 输出：新图（unique_ptr，所有权转移给调用方）
     * @param src  源图（只读）
     * @return SCIP_OKAY 成功，否则 SCIP_NOMEMORY
     */
    static SCIP_RETCODE copyGraph(std::unique_ptr<Graph> &dest, const Graph &src)
    {
        int num_nodes = src.numNodes();
        int num_edges = src.numEdges(); // 无向边数量

        // 创建新图：n 个节点，2*m 条半边（m 条无向边 × 2）
        dest = Graph::create(num_nodes, 2 * num_edges);
        if (!dest)
            return SCIP_NOMEMORY;

        // 复制节点：只复制静态属性，邻接表指针后续重建
        for (int i = 0; i < num_nodes; ++i)
        {
            auto &node = dest->node(i);
            const auto &srcnode = src.node(i);
            assert(srcnode.id == i);

            node.coord_x = srcnode.coord_x;
            node.coord_y = srcnode.coord_y;
            node.id = srcnode.id;
            node.first_edge = nullptr; // 邻接表将在下面重建
        }

        // 复制边：按原始遍历顺序重建邻接表和半边配对关系
        int edge_idx = 0;
        for (int i = 0; i < num_nodes - 1; ++i)
        {
            auto &start_node = dest->node(i);
            for (int j = i + 1; j < num_nodes; ++j)
            {
                auto &end_node = dest->node(j);

                // 半边对：k 为正向 (i→j)，k+num_edges 为反向 (j→i)
                auto &forward_edge = dest->edge(edge_idx);
                auto &reverse_edge = dest->edge(edge_idx + num_edges);

                // 设置目标节点
                forward_edge.target = &end_node;
                reverse_edge.target = &start_node;

                // 建立双向引用关系
                forward_edge.setReverseEdge(&reverse_edge);

                // 复制边长
                forward_edge.length = src.edge(edge_idx).length;
                reverse_edge.length = forward_edge.length;

                // 插入邻接表
                forward_edge.linkTo(&start_node);
                reverse_edge.linkTo(&end_node);

                ++edge_idx;
            }
        }

        return SCIP_OKAY;
    }

    /**
     * @brief 拷贝问题数据到另一个 SCIP 实例
     *
     * 深拷贝图数据，并将原始变量映射到目标 SCIP 实例的变量空间。
     * 用于 SCIP 的并行求解或子问题拷贝。
     */
    SCIP_RETCODE ProbDataTSP::scip_copy(
        SCIP *scip,
        SCIP *sourcescip,
        SCIP_HASHMAP *varmap,
        SCIP_HASHMAP *consmap,
        ObjProbData **objprobdata,
        SCIP_Bool global,
        SCIP_RESULT *result)
    {
        // 获取源数据
        ProbDataTSP *source_data = dynamic_cast<ProbDataTSP *>(SCIPgetObjProbData(sourcescip));
        assert(source_data != nullptr);
        const Graph *source_graph = source_data->getGraph();
        assert(source_graph != nullptr);

        // 深拷贝图结构
        std::unique_ptr<Graph> graph;
        SCIP_CALL(copyGraph(graph, *source_graph));

        // 将原始变量映射到目标 SCIP 实例中的对应变量
        int num_edges = graph->numEdges();
        for (int e = 0; e < num_edges; ++e)
        {
            auto &forward_edge = graph->edge(e);
            auto &reverse_edge = graph->edge(e + num_edges);
            assert(source_graph->edge(e).scip_var != nullptr);

            SCIP_Bool success;

            // 在目标 SCIP 实例中查找对应的变量
            SCIP_CALL(SCIPgetVarCopy(sourcescip, scip, source_graph->edge(e).scip_var,
                                     &forward_edge.scip_var, varmap, consmap, global, &success));
            // 增加变量引用计数（边持有该变量的引用）
            SCIP_CALL(SCIPcaptureVar(scip, forward_edge.scip_var));
            assert(success);
            assert(forward_edge.scip_var != nullptr);

            // 反向半边共享同一个变量
            reverse_edge.scip_var = forward_edge.scip_var;
            SCIP_CALL(SCIPcaptureVar(scip, reverse_edge.scip_var));
        }

        // 创建新的问题数据对象
        auto *probdata = new ProbDataTSP(std::move(graph));
        assert(probdata != nullptr);

        *objprobdata = probdata;
        *result = SCIP_SUCCESS;

        return SCIP_OKAY;
    }

    /**
     * @brief 释放原始问题数据
     *
     * 先释放所有边关联的 SCIP 变量（减少引用计数），再释放图本身。
     * releaseVar 会同时释放正向和反向半边的变量引用。
     */
    SCIP_RETCODE ProbDataTSP::scip_delorig(SCIP *scip)
    {
        for (int i = 0; i < graph_->numEdges(); ++i)
        {
            SCIP_CALL(graph_->edge(i).releaseVar(scip));
        }
        graph_.reset();
        return SCIP_OKAY;
    }

    /**
     * @brief 释放变换后的问题数据
     *
     * 与 scip_delorig 逻辑相同。在求解结束后释放变换数据。
     */
    SCIP_RETCODE ProbDataTSP::scip_deltrans(SCIP *scip)
    {
        for (int i = 0; i < graph_->numEdges(); ++i)
        {
            SCIP_CALL(graph_->edge(i).releaseVar(scip));
        }
        graph_.reset();
        return SCIP_OKAY;
    }

    /**
     * @brief 将原始问题数据变换为 presolved 问题数据
     *
     * 创建图的深拷贝，并将原始变量映射到变换后的变量。
     * 变换后的数据用于 presolve 后的求解阶段。
     */
    SCIP_RETCODE ProbDataTSP::scip_trans(
        SCIP *scip,
        ObjProbData **objprobdata,
        SCIP_Bool *deleteobject)
    {
        assert(objprobdata != nullptr);
        assert(deleteobject != nullptr);
        assert(graph_ != nullptr);

        // 深拷贝图结构
        std::unique_ptr<Graph> trans_graph;
        SCIP_CALL(copyGraph(trans_graph, *graph_));

        // 映射变量：将原始变量替换为变换后的对应变量
        int num_edges = trans_graph->numEdges();
        for (int e = 0; e < num_edges; ++e)
        {
            auto &forward_edge = trans_graph->edge(e);
            auto &reverse_edge = trans_graph->edge(e + num_edges);
            assert(graph_->edge(e).scip_var != nullptr);

            // 获取变换后的变量
            SCIP_CALL(SCIPgetTransformedVar(scip, graph_->edge(e).scip_var, &forward_edge.scip_var));
            // 增加引用计数
            SCIP_CALL(SCIPcaptureVar(scip, forward_edge.scip_var));

            // 反向半边共享变量
            reverse_edge.scip_var = forward_edge.scip_var;
            assert(reverse_edge.scip_var != nullptr);
            SCIP_CALL(SCIPcaptureVar(scip, reverse_edge.scip_var));
        }

        // 创建变换后的问题数据
        auto *trans_probdata = new ProbDataTSP(std::move(trans_graph));
        assert(trans_probdata != nullptr);

        *objprobdata = trans_probdata;
        *deleteobject = FALSE; // 保留原始数据：变换后的约束处理器通过 SCIP_ConsData
                               // 持有对原始图的引用，变换后仍需访问

        return SCIP_OKAY;
    }

}