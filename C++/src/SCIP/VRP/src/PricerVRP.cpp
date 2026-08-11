#include "PricerVRP.h"
#include "PQUEUE.h"

#include <iostream>
#include <map>
#include <vector>

#include "scip/cons_linear.h"

using namespace std;
using namespace scip;

/** Constructs the pricer object with the data needed
 *
 *  An alternative is to have a problem data class which allows to access the data.
 */
ObjPricerVRP::ObjPricerVRP(
    SCIP *scip,                                   /**< SCIP pointer */
    const char *p_name,                           /**< name of pricer */
    const int p_num_nodes,                        /**< number of nodes */
    const int p_capacity,                         /**< vehicle capacity */
    const vector<int> &p_demand,                  /**< demand array */
    const vector<vector<int>> &p_distance,        /**< matrix of distances */
    const vector<vector<SCIP_VAR *>> &p_arc_var,  /**< matrix of arc variables */
    const vector<vector<SCIP_CONS *>> &p_arc_con, /**< matrix of arc constraints */
    const vector<SCIP_CONS *> &p_part_con         /**< array of partitioning constraints */
    ) : ObjPricer(scip, p_name, "Finds tour with negative reduced cost.", 0, TRUE),
        _num_nodes(p_num_nodes),
        _capacity(p_capacity),
        _demand(p_demand),
        _distance(p_distance),
        _arc_var(p_arc_var),
        _arc_con(p_arc_con),
        _part_con(p_part_con)
{
}

/* Destructs the pricer object */
ObjPricerVRP::~ObjPricerVRP() {}

/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerVRP::scip_init)
{
   for (int i = 0; i < num_nodes(); ++i)
   {
      if (i > 0)
      {
         SCIP_CALL(SCIPgetTransformedCons(scip, _part_con[i], &_part_con[i]));
      }

      for (int j = 0; j < i; ++j)
      {
         SCIP_CALL(SCIPgetTransformedVar(scip, _arc_var[i][j], &_arc_var[i][j]));
         SCIP_CALL(SCIPgetTransformedCons(scip, _arc_con[i][j], &_arc_con[i][j]));
      }
   }

   return SCIP_OKAY;
}

/**
 * perform pricing
 *
 * @todo compute shortest length restricted tour w.r.t duals
 */
SCIP_RETCODE ObjPricerVRP::pricing(SCIP *scip, bool isfarkas) const
{
   int numNodes = num_nodes();
   /* 创建检验数数组 */
   vector<vector<SCIP_Real>> red_length(numNodes);
   for (int i = 0; i < numNodes; ++i)
   {
      red_length[i].resize(i, 0.0);
   }

   /* 计算缩减成本的弧长，只存储下三角矩阵，即仅当 i > j 时存储 red_length[i][j] */
   if (isfarkas)
   {
      for (int i = 0; i < numNodes; ++i)
      {
         assert(i == 0 || part_con(i) != 0);
         for (int j = 0; j < i; ++j)
         {
            SCIP_Real r = 0.0;
            assert(arc_con(i, j) != 0);

            r -= SCIPgetDualfarkasLinear(scip, arc_con(i, j));
            if (j != 0)
            {
               r -= 0.5 * SCIPgetDualfarkasLinear(scip, part_con(j));
            }
            if (i != 0)
            {
               r -= 0.5 * SCIPgetDualfarkasLinear(scip, part_con(i));
            }
            red_length[i][j] = r;
         }
      }
   }
   else
   {
      for (int i = 0; i < numNodes; ++i)
      {
         assert(i == 0 || part_con(i) != 0);
         for (int j = 0; j < i; ++j)
         {
            SCIP_Real r = 0.0;
            assert(arc_con(i, j) != 0);
            r -= SCIPgetDualsolLinear(scip, arc_con(i, j));

            if (j != 0)
               r -= 0.5 * SCIPgetDualsolLinear(scip, part_con(j));
            if (i != 0)
               r -= 0.5 * SCIPgetDualsolLinear(scip, part_con(i));
            red_length[i][j] = r;
         }
      }
   }
#ifdef SCIP_OUTPUT
   if (isfarkas)
   {
      SCIPinfoMessage(scip, NULL, "dual ray solution:\n");
      for (int i = 0; i < num_nodes(); ++i)
      {
         for (int j = 0; j < i; ++j)
            SCIPinfoMessage(scip, NULL, "arc_%d_%d:  %g\n", i, j, SCIPgetDualfarkasLinear(scip, arc_con(i, j)));
      }

      for (int i = 1; i < num_nodes(); ++i)
         SCIPinfoMessage(scip, NULL, "part_%d:  %g\n", i, SCIPgetDualfarkasLinear(scip, part_con(i)));

      for (int i = 0; i < num_nodes(); ++i)
      {
         for (int j = 0; j < i; ++j)
            SCIPinfoMessage(scip, NULL, "length_%d_%d:  %g\n", i, j, red_length[i][j]);
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "dual solution:\n");
      for (int i = 0; i < num_nodes(); ++i)
      {
         for (int j = 0; j < i; ++j)
            SCIPinfoMessage(scip, NULL, "arc_%d_%d:  %g\n", i, j, SCIPgetDualsolLinear(scip, arc_con(i, j)));
      }

      for (int i = 1; i < num_nodes(); ++i)
         SCIPinfoMessage(scip, NULL, "part_%d:  %g\n", i, SCIPgetDualsolLinear(scip, part_con(i)));

      for (int i = 0; i < num_nodes(); ++i)
      {
         for (int j = 0; j < i; ++j)
            SCIPinfoMessage(scip, NULL, "length_%d_%d:  %g\n", i, j, red_length[i][j]);
      }
   }
#endif

   /* 使用动态规划计算检验数 */
   list<int> tour;
   SCIP_Real reduced_cost = find_shortest_tour(red_length, tour);

   /* 向主问题中添加变量 */
   if (SCIPisNegative(scip, reduced_cost))
   {
      return add_tour_variable(scip, tour);
   }

#ifdef SCIP_OUTPUT
   SCIP_CALL(SCIPwriteTransProblem(scip, "vrp.lp", "lp", FALSE));
#endif

   return SCIP_OKAY;
}

/** 当 LP 可行时，对附加变量进行定价。
 *
 *  - 获取所需对偶变量的值
 *  - 根据这些值构建具有缩减成本的弧长（reduced-cost arc lengths）
 *  - 寻找在这些弧长下最短的合法回路（admissible tour）
 *  - 如果该回路的缩减成本为负，则将其添加到 LP 中
 *
 *  *result 的可能返回值：
 *  - SCIP_SUCCESS    : 至少找到了一个能改善目标的变量，或者已确保不存在这样的变量
 *  - SCIP_DIDNOTRUN  : 定价过程被定价器中止，无法保证当前的 LP 解是最优的
 */
SCIP_DECL_PRICERREDCOST(ObjPricerVRP::scip_redcost)
{
   SCIPinfoMessage(scip, NULL, "call scip_redcost ...\n");

   /* set result pointer, see above */
   *result = SCIP_SUCCESS;

   /* call pricing routine */
   SCIP_CALL(pricing(scip, false));

   return SCIP_OKAY;
}

/** 当 LP 不可行时，对附加变量进行定价。
 *
 *  - 获取所需的对偶 Farkas 乘子（dual Farkas multipliers）的值
 *  - 根据这些值构建具有缩减成本的弧长
 *  - 寻找在这些弧长下最短的合法回路
 *  - 如果该回路的缩减成本为负，则将其添加到 LP 中
 */
SCIP_DECL_PRICERFARKAS(ObjPricerVRP::scip_farkas)
{
   SCIPinfoMessage(scip, NULL, "call scip_farkas ...\n");

   /* call pricing routine */
   SCIP_CALL(pricing(scip, true));

   return SCIP_OKAY;
} /*lint !e715*/

/* 向主问题中添加变量 */
SCIP_RETCODE ObjPricerVRP::add_tour_variable(SCIP *scip, const list<int> &tour) const
{
   char tmp_name[255];
   char var_name[255];
   (void)SCIPsnprintf(var_name, 255, "T");
   for (list<int>::const_iterator it = tour.begin(); it != tour.end(); ++it)
   {
      (void)strncpy(tmp_name, var_name, 255);
      (void)SCIPsnprintf(var_name, 255, "%s_%d", tmp_name, *it);
   }
   SCIPinfoMessage(scip, NULL, "new variable <%s>\n", var_name);

   /* 创建新变量：将上界设为无穷大，这样在定价过程中就无需关心该变量的缩减成本了。
    * 上界为 1 的约束会通过集合划分约束（set partitioning constraints）隐式地得到满足。
    */
   SCIP_VAR *var;
   SCIP_CALL(SCIPcreateVar(scip, &var, var_name,
                           0.0,
                           SCIPinfinity(scip),
                           0.0,
                           SCIP_VARTYPE_CONTINUOUS,
                           FALSE, FALSE, NULL, NULL, NULL, NULL, NULL));

   SCIP_CALL(SCIPaddPricedVar(scip, var, 1.0));

   /* 向集合划分约束中添加变量以及系数 */
   for (list<int>::const_iterator it = tour.begin(); it != tour.end(); ++it)
   {
      assert(0 <= *it && *it <= num_nodes());
      SCIP_CALL(SCIPaddCoefLinear(scip, part_con(*it), var, 1.0));
   }

   /* 向弧-路径耦合约束中添加变量 */
   int last = 0;
   for (list<int>::const_iterator it = tour.begin(); it != tour.end(); ++it)
   {
      assert(0 <= *it && *it < num_nodes());
      SCIP_CALL(SCIPaddCoefLinear(scip, arc_con(last, *it), var, 1.0));
      last = *it;
   }
   SCIP_CALL(SCIPaddCoefLinear(scip, arc_con(last, 0), var, 1.0));

   /* cleanup */
   SCIP_CALL(SCIPreleaseVar(scip, &var));

   return SCIP_OKAY;
}

/** 根据给定的距离（或长度）计算出一条最短的可行路线。
 *  函数必须通过 tour 参数返回计算出的路线，并将该路线的长度（基于给定的距离）作为返回值。
 *  返回的路线必须是一个包含路线中所有客户节点的有序列表
 *  （例如：对于 0-2-5-7-0 这样一条路线，应返回 2-5-7）。
 */

namespace
{
   constexpr SCIP_Real eps = 1e-9;

   struct PQUEUE_KEY
   {
      int demand;
      SCIP_Real length;

      PQUEUE_KEY() : demand(0), length(0.0) {}
   };

   bool operator<(const PQUEUE_KEY &l1, const PQUEUE_KEY &l2)
   {
      if (l1.demand < l2.demand)
      {
         return true;
      }
      if (l1.demand > l2.demand)
      {
         return false;
      }
      if (l1.length < l2.length - eps)
      {
         return true;
      }
      return false;
   }

   using PQUEUE_DATA = int;
   using PQUEUE = pqueue<PQUEUE_KEY, PQUEUE_DATA>;
   using PQUEUE_ITEM = PQUEUE::pqueue_item;

   struct NODE_TABLE_DATA
   {
      SCIP_Real length;
      int prodecessor;
      PQUEUE_ITEM queue_item;

      NODE_TABLE_DATA() : length(0.0), prodecessor(-1), queue_item(nullptr) {}
   };

   using NODE_TABLE_KEY = int;
   using NODE_TABLE = std::map<NODE_TABLE_KEY, NODE_TABLE_DATA>;

}

/** 返回负约化成本回路（使用受限最短路径动态规划算法）
 *
 *  该算法使用了 pqueue.h 中的优先队列实现。之所以不能用 SCIP 自带的优先队列，
 *  是因为它目前还不支持删除非队首（非顶部）的元素。
 */
SCIP_Real ObjPricerVRP::find_shortest_tour(const vector<vector<SCIP_Real>> &length, list<int> &tour) const
{
   tour.clear();

   SCIPdebugMessage("Enter RSP - capacity: %d\n", capacity());

   /* 算法开始 */
   PQUEUE PQ;
   
   vector<NODE_TABLE> table(num_nodes());

   /* 插入根节点 */
   PQUEUE_KEY queue_key;
   PQUEUE_DATA queue_data = 0;
   PQUEUE_ITEM queue_item = PQ.insert(queue_key, queue_data);

   NODE_TABLE_KEY table_key = 0;
   NODE_TABLE_DATA table_entry;

   while (!PQ.empty())
   {
      queue_item = PQ.top();
      queue_key = PQ.get_key(queue_item);
      queue_data = PQ.get_data(queue_item);
      PQ.pop();

      const int curr_node = queue_data;
      const SCIP_Real curr_length = queue_key.length;
      const int curr_demand = queue_key.demand;

      /* 节点为根节点，并且路径为负值，说明找到了可行路径*/
      if (curr_node == 0 && curr_length < -eps)
      {
         break;
      }

      if (curr_node == 0 && curr_demand != 0)
      {
         continue;
      }

      for (int next_node = 0; next_node < num_nodes(); ++next_node)
      {
         if (next_node == curr_node)
         {
            continue;
         }

         if (have_edge(next_node, curr_node) == false)
         {
            continue;
         }

         const int next_demand = curr_demand + demand(next_node);

         if (next_demand > capacity())
         {
            continue;
         }

         const SCIP_Real next_length = curr_length + (curr_node > next_node ? length[curr_node][next_node] : length[next_node][curr_node]);
         NODE_TABLE &next_table = table[next_node];

         bool skip = false;
         list<NODE_TABLE::iterator> dominated;

         for (auto it = next_table.begin(); it != next_table.end() && !skip; ++it)
         {
            if (next_demand >= it->first && next_length >= it->second.length - eps)
            {
               skip = true;
            }

            if (next_demand <= it->first && next_length <= it->second.length - eps)
            {
               dominated.push_front(it);
            }
         }
         if (skip)
         {
            continue;
         }

         for (auto it = dominated.begin(); it != dominated.end(); ++it)
         {
            PQ.remove((*it)->second.queue_item);
            (void)next_table.erase(*it);
         }

         /* 插入新的节点，更新table*/
         queue_key.demand = next_demand;
         queue_key.length = next_length;
         queue_data = next_node;

         queue_item = PQ.insert(queue_key, queue_data);

         table_key = next_demand;
         table_entry.length = next_length;
         table_entry.prodecessor = curr_node;
         table_entry.queue_item = queue_item;

         next_table[table_key] = table_entry;
#ifdef SCIP_OUTPUT
         printf("new entry  node = %d  demand = %d  length = %g  pref = %d\n", next_node, next_demand, next_length, curr_node);
#endif
      }
   }

   SCIPdebugMessage("Done RSP DP.\n");

   table_entry.prodecessor = -1;
   table_entry.length = 0;
   int curr_node = 0;

   for (auto it = table[0].begin(); it != table[0].end(); ++it)
   {
      if (it->second.length < table_entry.length) 
      {
         table_key = it->first;
         table_entry = it->second;
      }
   }
   SCIP_Real tour_length = table_entry.length;
   
   while (table_entry.prodecessor > 0)
   {
      table_key -= demand(curr_node);
      curr_node = table_entry.prodecessor;
      tour.push_front(curr_node);
      table_entry = table[curr_node][table_key];
   }
   SCIPdebugMessage("Leave RSP  tour length = %g\n", tour_length);
   return tour_length;
}