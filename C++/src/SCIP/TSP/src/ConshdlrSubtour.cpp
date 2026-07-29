/**@file   ConshdlrSubtour.cpp
 * @brief  TSP 子环消除约束处理器实现
 * @author Timo Berthold
 *
 * @details
 * 实现 TSP 子环消除约束的所有 SCIP 回调函数。
 *
 * ## 核心流程
 *
 * ### 分离（Separation）
 * 1. 从 LP 解中读取每条边的变量值，设置为图的边容量
 * 2. 调用 ghc_tree() 构建 Gomory-Hu 割树，获取所有违反的割
 * 3. 将每个割转换为 SCIP 行约束（Σ_{u∈S, v∉S} x_uv ≥ 2）
 * 4. 检查割是否有效（对当前 LP 解有违反），有效则添加到 SCIP
 *
 * ### 可行性检查（Check）
 * 从节点 0 出发沿 LP 解中选中的边遍历图。如果能回到节点 0
 * 且访问节点数 = 总节点数，则解可行（完整哈密顿环）；否则存在子环。
 *
 * ## 子环检测原理
 *
 * 在 LP 松弛解中，变量 x_e 取值在 [0,1] 之间。即使度约束满足
 * （每个节点出入度和 = 2），LP 解仍可能表示多个不相交的子环
 * （子环 / subtour）。例如：
 *
 * ~~~
 *   6 个节点的图中，LP 解可能是：
 *   子环 1: 0→1→2→0  (x 值 = 1)
 *   子环 2: 3→4→5→3  (x 值 = 1)
 *   总度数 = 2，但没有形成单一哈密顿环
 * ~~~
 *
 * 通过 Gomory-Hu 树可以找到"连接较弱"的节点子集，
 * 生成割平面来切断这些子环。
 */

#include <cassert>
#include <string>
#include <iostream>
#include "ConshdlrSubtour.h"
#include "GomoryHuTree.h"

#include "objscip/objscip.h"
#include "scip/cons_linear.h"

using namespace scip;

namespace tsp
{

   /* ==========================================================================
    * 辅助函数
    * ========================================================================== */

   /**
    * @brief 将当前 LP 解的变量值同步为图中边的容量
    *
    * 在每次调用 ghc_tree() 之前执行。
    * Gomory-Hu 树算法需要边容量作为最大流计算的输入。
    * 这里用 LP 解的当前变量值作为容量（每条边对应一个二元变量 x_e ∈ [0,1]，
    * LP 松弛后为小数值）。
    *
    * 边容量决定了在 Gomory-Hu 树中哪些部分"连接紧密"（容量大）
    * 哪些部分"连接较弱"（容量小），从而找到可以生成割平面的位置。
    *
    * @param scip   SCIP 实例
    * @param sol    当前 LP 解
    * @param graph  图（边的 capacity 和 residual_capacity 将被设置为 LP 值）
    */
   static void resetEdgeCapacitiesFromLP(SCIP *scip, SCIP_SOL *sol, Graph *graph)
   {
      int num_edges = graph->numEdges();

      for (int i = 0; i < num_edges; ++i)
      {
         auto &edge = graph->edge(i);
         assert(edge.scip_var != nullptr);
         double cap = SCIPgetSolVal(scip, sol, edge.scip_var);

         // 无向图：正向半边和反向半边对称设置
         edge.capacity = cap;
         edge.residual_capacity = cap;
         edge.reverse->capacity = cap;
         edge.reverse->residual_capacity = cap;
      }
   }

   /**
    * @brief 从一个割（位向量）构造 SCIP 行约束
    *
    * 割 cuts[i] 是一个长度为 num_nodes 的位向量：
    *   - cuts[i][j] = true  → 节点 j 在汇侧 (T)
    *   - cuts[i][j] = false → 节点 j 在源侧 (S)
    *
    * 子环消除不等式的形式：
    * ~~~
    *   Σ_{u∈S, v∈T} x_{uv} ≥ 2
    * ~~~
    *
    * 即：跨越割的边数 >= 2（如果小于 2，说明子集 S 形成了隔离的子环）。
    *
    * @param scip      SCIP 实例
    * @param conshdlr  约束处理器
    * @param cut       割的位向量
    * @param graph     图
    * @param row       输出：构造的 SCIP 行（调用方负责释放）
    * @return SCIP_OKAY 成功，否则错误码
    */
   static SCIP_RETCODE createRowFromCut(
       SCIP *scip,
       SCIP_CONSHDLR *conshdlr,
       const std::vector<bool> &cut,
       Graph *graph,
       SCIP_ROW **row)
   {
      int num_nodes = graph->numNodes();

      SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, row, conshdlr, "subtour",
                                           2.0, SCIPinfinity(scip), FALSE, FALSE, TRUE));
      SCIP_CALL(SCIPcacheRowExtensions(scip, *row));

      for (int j = 0; j < num_nodes; ++j)
      {
         if (!cut[j])
         {
            continue;
         }
         auto &node = graph->node(j);
         for (auto *edge = node.first_edge; edge != nullptr; edge = edge->next)
         {
            if (!cut[edge->target->id])
            {
               SCIP_CALL(SCIPaddVarToRow(scip, *row, edge->scip_var, 1.0));
            }
         }
      }

      SCIP_CALL(SCIPflushRowExtensions(scip, *row));
      return SCIP_OKAY;
   }

   /**
    * @brief 在节点的邻接表中查找下一条有效的遍历边
    *
    * 在 TSP 解的图中，从当前节点出发，查找一条被 LP 解选中
    * （变量值 > 0.5）且不是刚刚走过的边（跳过反向边）的邻接边。
    *
    * @param node       当前节点
    * @param last_edge  上一步走过的边，其反向边将被跳过
    * @param scip       SCIP 实例
    * @param sol        当前解
    * @return 找到的有效边指针；若没有则返回 nullptr（路径中断，存在子环）
    */
   static GraphEdge *findValidNextEdge(
       GraphNode *node,
       GraphEdge *last_edge,
       SCIP *scip,
       SCIP_SOL *sol)
   {
      for (auto *edge = node->first_edge; edge != nullptr; edge = edge->next)
      {
         // 跳过刚走过的边（如果立即回头，说明没有形成简单路径）
         if (edge->reverse == last_edge)
            continue;

         // 检查解中该边是否被选中：LP 解变量值是连续的，
         // 整数解变量值是二值的，> 0.5 判定为选中
         if (SCIPgetSolVal(scip, sol, edge->scip_var) > 0.5)
            return edge;
      }
      return nullptr;
   }
   /**
    * @brief 检查当前解是否包含子环
    *
    * 从节点 0 出发，沿解中选中的边遍历图：
    *   - 如果能回到节点 0 且访问节点数 == 总节点数 → 完整哈密顿环，可行
    *   - 如果路径中断或访问节点数 < 总节点数 → 存在子环，不可行
    *
    * @param scip   SCIP 实例
    * @param graph  TSP 图
    * @param sol    当前解
    * @return true 如果发现子环（解不可行），false 如果形成完整环（解可行）
    */
   static bool findSubtour(SCIP *scip, Graph *graph, SCIP_SOL *sol)
   {
      assert(scip != nullptr);
      assert(graph != nullptr);

      int num_nodes = graph->numNodes();
      if (num_nodes <= 1)
      {
         return false;
      }

      auto *current_node = &graph->node(0);
      auto *start_node = current_node;
      GraphEdge *last_edge = nullptr;
      int visited_count = 0;

      do
      {
         auto *next_edge = findValidNextEdge(current_node, last_edge, scip, sol);

         // 路径中断 -> 存在子环
         if (next_edge == nullptr || visited_count > num_nodes)
         {
            return true;
         }
         ++visited_count;
         last_edge = next_edge;
         current_node = next_edge->target;
      } while (current_node != start_node);
      return visited_count != num_nodes;
   }
   /* ==========================================================================
    * 分离回调
    *
    * 当 LP 松弛解被求解后，SCIP 调用此回调来寻找该解违反的约束。
    * 如果找到，添加到 LP 中以切断当前非法解。
    * ========================================================================== */

   /**
    * @brief LP 分离 / 整数解分离 / 强制执行的回调
    *
    * 对每个有效约束调用 Ghomory-Hu 树算法生成割平面。
    *
    * @param is_enforcement 是否强制模式（true = 必须找到有效割或返回 CUTOFF）
    */
   static SCIP_RETCODE sepaSubtour(
       SCIP *scip,
       SCIP_CONSHDLR *conshdlr,
       SCIP_CONS **constraints,
       int num_constraints,
       int num_useful_constraints,
       SCIP_SOL *sol,
       bool is_enforcement,
       SCIP_RESULT *result)
   {
      assert(result != nullptr);
      *result = SCIP_DIDNOTFIND;

      // 对每个活跃约束执行分离
      for (int c = 0; c < num_useful_constraints && *result != SCIP_CUTOFF; ++c)
      {
         SCIP_ConsData *constraint_data = SCIPconsGetData(constraints[c]);
         assert(constraint_data != nullptr);

         Graph *graph = constraint_data->graph;
         assert(graph != nullptr);

         // 步骤1：将LP解的变量值设置为边的容量
         resetEdgeCapacitiesFromLP(scip, sol, graph);

         // 步骤2：构建 Gomory_Hu 割树并获取所有违反的割
         std::vector<std::vector<bool>> cuts;
         int num_cuts;

         if (ghc_tree(graph, cuts, &num_cuts, SCIPfeastol(scip)))
         {
            // 步骤3：对每个违反的割，构造SCIP约束
            for (int i = 0; i < num_cuts && *result != SCIP_CUTOFF; i++)
            {
               SCIP_ROW *row;
               SCIP_CALL(createRowFromCut(scip, conshdlr, cuts[i], graph, &row));

               // 步骤4：判断割是否有效
               // 在强制模式下直接添加，否则需要检查割的“效用”->对LP解的违反程度
               if (is_enforcement || SCIPisCutEfficacious(scip, sol, row))
               {
                  SCIP_Bool infeasible;
                  SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));

                  // 如果割不可行(CUT_OFF)，问题在该分支无解
                  *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
               }
               SCIP_CALL(SCIPreleaseRow(scip, &row));
            }
         }
      }

      return SCIP_OKAY;
   }

   /* ==========================================================================
    * SCIP 回调实现
    *
    * 以下是 ConshdlrSubtour 类的所有虚函数实现。
    * 每个函数对应 SCIP 约束处理器框架的一个回调点。
    * ========================================================================== */
   /**
    * @brief 约束删除回调
    *
    * 约束被 SCIP 释放时调用。图数据的所有权归 ProbDataTSP，
    * 这里不需要释放 graph 指针。
    */
   SCIP_DECL_CONSDELETE(ConshdlrSubtour::scip_delete)
   { /*lint --e{715}*/
      assert(scip != nullptr);
      assert(cons != nullptr);
      assert(consdata != nullptr);
      assert(*consdata != nullptr);

      // Graph 的所有权在 ProbDataTSP 中，这里只需标记删除
      delete *consdata;
      *consdata = nullptr;

      return SCIP_OKAY;
   }

   /**
    * @brief 约束变换回调
    *
    * 问题从原始空间变换到 presolved 空间时调用。
    * 需要深拷贝约束数据（特别是图结构），因为 SCIP 在变换后
    * 会释放原始问题数据。
    */
   SCIP_DECL_CONSTRANS(ConshdlrSubtour::scip_trans)
   { /*lint --e{715}*/
      assert(scip != nullptr);
      assert(sourcecons != nullptr);
      assert(targetcons != nullptr);

      SCIP_ConsData *source_data = SCIPconsGetData(sourcecons);
      assert(source_data != nullptr);

      // 与 C 版本一致：创建新的约束数据，指向原始图
      // SCIP 内部会自动将原始变量映射为变换后的变量
      // 注意：原始图必须保持存活，不能随 ProbDataTSP::scip_delorig 释放
      auto *target_data = new SCIP_ConsData(source_data->graph);

      SCIP_CALL(SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons),
                               conshdlr, target_data,
                               SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
                               SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
                               SCIPconsIsPropagated(sourcecons),
                               SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
                               SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
                               SCIPconsIsStickingAtNode(sourcecons)));

      return SCIP_OKAY;
   }

   /**
    * @brief LP 分离回调
    *
    * 在 LP 松弛求解后被调用。使用 Gomory-Hu 树查找当前 LP 解违反的子环消除不等式。
    */
   SCIP_DECL_CONSSEPALP(ConshdlrSubtour::scip_sepalp)
   { /*lint --e{715}*/
      *result = SCIP_DIDNOTFIND;

      // 调用通用分离函数，enforce = false（不强制）
      SCIP_CALL(sepaSubtour(scip, conshdlr, conss, nconss, nusefulconss,
                            nullptr, false, result));

      return SCIP_OKAY;
   }

   /**
    * @brief 整数解分离回调
    *
    * 当找到整数可行解时调用。沿图遍历检查是否构成完整哈密顿环。
    * 如果发现子环，返回 DIDNOTFIND 让 enforce 回调处理。
    */
   SCIP_DECL_CONSSEPASOL(ConshdlrSubtour::scip_sepasol)
   {
      *result = SCIP_FEASIBLE;

      for (int c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c)
      {
         SCIP_ConsData *constraint_data = SCIPconsGetData(conss[c]);
         assert(constraint_data != nullptr);

         // 检查解中是否存在子环
         if (findSubtour(scip, constraint_data->graph, sol))
         {
            *result = SCIP_DIDNOTFIND;
         }
      }
      return SCIP_OKAY;
   }

   /**
    * @brief LP 强制执行回调
    *
    * 当 LP 解违反约束且需要强制切断时调用。
    * 与 sepalp 类似但 enforce = true（必须找到有效割或返回不可行）。
    */
   SCIP_DECL_CONSENFOLP(ConshdlrSubtour::scip_enfolp)
   { /*lint --e{715}*/
      *result = SCIP_FEASIBLE;

      // 调用分离函数，enforce = true（强制模式）
      SCIP_CALL(sepaSubtour(scip, conshdlr, conss, nconss, nusefulconss,
                            nullptr, true, result));

      return SCIP_OKAY;
   }

   /**
    * @brief 伪解强制执行回调
    *
    * 当 SCIP 的伪解（pseudo solution，所有整数变量置为边界值）违反约束时调用。
    * 对 TSP 子环约束，检查伪解是否构成完整哈密顿环。
    */
   SCIP_DECL_CONSENFOPS(ConshdlrSubtour::scip_enfops)
   { /*lint --e{715}*/
      *result = SCIP_FEASIBLE;

      for (int c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c)
      {
         SCIP_ConsData *constraint_data = SCIPconsGetData(conss[c]);
         assert(constraint_data != nullptr);

         if (findSubtour(scip, constraint_data->graph, nullptr))
            *result = SCIP_DIDNOTFIND;
      }

      return SCIP_OKAY;
   }

   /**
    * @brief 可行性检查回调
    *
    * 验证给定解是否满足子环消除约束。
    * 检查解是否构成一个完整的哈密顿环（而非多个子环）。
    */
   SCIP_DECL_CONSCHECK(ConshdlrSubtour::scip_check)
   { /*lint --e{715}*/
      *result = SCIP_FEASIBLE;

      for (int c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c)
      {
         SCIP_ConsData *constraint_data = SCIPconsGetData(conss[c]);
         assert(constraint_data != nullptr);

         if (findSubtour(scip, constraint_data->graph, sol))
            *result = SCIP_INFEASIBLE;
      }

      return SCIP_OKAY;
   }

   /**
    * @brief 约束传播回调
    *
    * 在 presolve 和 constraint propagation 期间被调用。
    * 对每个约束处理传播逻辑。TSP 子环约束的传播能力有限，
    * 主要依赖分离阶段生成割平面。
    */
   SCIP_DECL_CONSPROP(ConshdlrSubtour::scip_prop)
   { /*lint --e{715}*/
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   /**
    * @brief 变量锁定回调
    *
    * 锁定约束涉及的所有变量，防止它们在 rounding 时被修改。
    * 对于子环约束，锁定图中所有边的变量。
    */
   SCIP_DECL_CONSLOCK(ConshdlrSubtour::scip_lock)
   { /*lint --e{715}*/
      // SCIP_DECL_CONSLOCK 参数：scip, cons, consdata, locktype, nlockspos, nlocksneg
      SCIP_ConsData *constraint_data = SCIPconsGetData(cons);
      assert(constraint_data != nullptr);

      Graph *graph = constraint_data->graph;
      int num_edges = graph->numEdges();

      // 锁定所有边的变量（每条无向边对应一个 SCIP 变量）
      for (int i = 0; i < num_edges; ++i)
      {
         auto &edge = graph->edge(i);
         SCIP_CALL(SCIPaddVarLocksType(scip, edge.scip_var, locktype, nlockspos, nlocksneg));
      }

      return SCIP_OKAY;
   }

   /**
    * @brief 变量删除回调
    *
    * 当约束涉及的变量被删除时调用。TSP 子环约束的边变量不应被删除，
    * 此回调为占位。
    */
   SCIP_DECL_CONSDELVARS(ConshdlrSubtour::scip_delvars)
   { /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /**
    * @brief 约束打印回调
    *
    * 将约束输出到文件（用于调试和日志）。
    * 打印图中的边和对应的变量名。
    */
   SCIP_DECL_CONSPRINT(ConshdlrSubtour::scip_print)
   { /*lint --e{715}*/
      // SCIP_DECL_CONSPRINT 参数：scip, conshdlr, cons, file
      // 注意：consdata 不是直接参数，需要通过 SCIPconsGetData() 获取
      assert(scip != nullptr);
      assert(cons != nullptr);

      SCIP_ConsData *constraint_data = SCIPconsGetData(cons);
      assert(constraint_data != nullptr);

      Graph *graph = constraint_data->graph;
      int num_edges = graph->numEdges();

      for (int i = 0; i < num_edges; ++i)
      {
         auto &edge = graph->edge(i);
         if (i > 0)
            SCIPinfoMessage(scip, file, ", ");

         // 打印变量名
         SCIP_CALL(SCIPwriteVarName(scip, file, edge.scip_var, TRUE));
      }

      return SCIP_OKAY;
   }

   /**
    * @brief 约束处理器克隆回调
    *
    * 当 SCIP 需要克隆约束处理器时调用（例如并行求解中的子 SCIP）。
    */
   SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable *ConshdlrSubtour::clone) /*lint !e665*/
   {                                                                 /*lint --e{715}*/
      // SCIP_DECL_CONSHDLRCLONE 宏展开为：ObjProbCloneable* ConshdlrSubtour::clone(SCIP* scip, SCIP_Bool* valid) const
      // 约束处理器是无状态的，标记为有效并返回新实例
      assert(valid != NULL);
      *valid = TRUE;
      return new ConshdlrSubtour(scip);
   }

   /**
    * @brief 约束拷贝回调
    *
    * 将约束拷贝到另一个 SCIP 实例中。
    */
   SCIP_DECL_CONSCOPY(ConshdlrSubtour::scip_copy)
   { /*lint --e{715}*/
      assert(scip != nullptr);
      assert(cons != nullptr);

      // 从目标 SCIP 实例获取问题数据中的图（ProbDataTSP::scip_copy 已深拷贝）
      ProbDataTSP *probdata = dynamic_cast<ProbDataTSP *>(SCIPgetObjProbData(scip));
      assert(probdata != nullptr);
      Graph *graph = probdata->getGraph();

      // 创建约束的副本（使用目标实例的图，而非 nullptr）
      SCIP_CALL(SCIPcreateConsSubtour(scip, cons, SCIPconsGetName(sourcecons),
                                      graph, initial, separate, enforce, check, propagate,
                                      local, modifiable, dynamic, removable));

      return SCIP_OKAY;
   }

   /* ==========================================================================
    * 约束创建函数
    * ========================================================================== */

   /**
    * @brief 创建并注册一个 TSP 子环消除约束
    *
    * 分配 SCIP_ConsData，设置 SCIP 约束属性，将约束添加到 SCIP 实例。
    */
   SCIP_RETCODE SCIPcreateConsSubtour(
       SCIP *scip,
       SCIP_CONS **cons,
       const char *name,
       Graph *graph,
       SCIP_Bool initial,
       SCIP_Bool separate,
       SCIP_Bool enforce,
       SCIP_Bool check,
       SCIP_Bool propagate,
       SCIP_Bool local,
       SCIP_Bool modifiable,
       SCIP_Bool dynamic,
       SCIP_Bool removable)
   {
      auto *consdata = new SCIP_ConsData(graph);

      SCIP_CALL(SCIPcreateCons(scip, cons, name,
                               SCIPfindConshdlr(scip, "subtour"), consdata,
                               initial, separate, enforce, check, propagate,
                               local, modifiable, dynamic, removable, FALSE));

      return SCIP_OKAY;
   }

} // end tsp