/**@file   ProbDataTSP.cpp
 * @brief  TSP 问题的 C++ 问题数据类实现
 * @author GaiGuoNing
 */

#include "objscip/objscip.h"
#include "ProbData.h"
#include "GomoryHuTree.h"

using namespace tsp;
using namespace scip;

/** 复制给定的图 */
static
SCIP_RETCODE copy_graph(
   GRAPH**               graph,              /**< 用于存储复制后图的指针 */
   GRAPH*                sourcegraph         /**< 要被复制的源图 */
   )
{
   assert( graph != NULL );
   assert( sourcegraph != NULL );

   // 按照文件读取器中创建图的方式复制图
   int n = sourcegraph->nnodes;
   int m = sourcegraph->nedges;

   // create_graphs 为每条边分配两个反向平行弧（anti-parallel arcs）的内存
   if( ! create_graph(n, 2*m, graph))
      return SCIP_NOMEMORY;

   // 复制节点
   for(int i = 0; i < n; ++i)
   {
      GRAPHNODE* node       = &((*graph)->nodes[i]);
      GRAPHNODE* sourcenode = &(sourcegraph->nodes[i]);

      assert(sourcenode->id == i);

      node->x          = sourcenode->x;
      node->y          = sourcenode->y;
      node->id         = sourcenode->id;
      node->first_edge = NULL;
   }

   // 复制边
   int e = 0;
   for(int i = 0; i < n - 1; ++i)
   {
      GRAPHNODE* nodestart = &((*graph)->nodes[i]);
      for(int j = i + 1; j < n; ++j)
      {
         GRAPHNODE* nodeend = &((*graph)->nodes[j]);
         GRAPHEDGE* edgeforw  = &((*graph)->edges[e]);
         GRAPHEDGE* edgebackw = &((*graph)->edges[e + m]);

         // 构造两个“平行”的半边缘（halfedges）
         edgeforw->adjac  = nodeend;
         edgebackw->adjac = nodestart;
         edgeforw->back   = edgebackw;
         edgebackw->back  = edgeforw;

         // 复制长度
         edgeforw->length  = sourcegraph->edges[e].length;
         edgebackw->length = edgeforw->length;

         // 将其中一个半边缘插入到起始节点的边列表中
         if (nodestart->first_edge == NULL)
         {
            nodestart->first_edge = edgeforw;
            nodestart->first_edge->next = NULL;
         }
         else
         {
            edgeforw->next = nodestart->first_edge;
            nodestart->first_edge = edgeforw;
         }

         // 同理，将另一个半边缘插入到结束节点的边列表中
         if (nodeend->first_edge == NULL)
         {
            nodeend->first_edge = edgebackw;
            nodeend->first_edge->next = NULL;
         }
         else
         {
            edgebackw->next = nodeend->first_edge;
            nodeend->first_edge = edgebackw;
         }

         ++e;
      } // for j
   } // for i

   return SCIP_OKAY;
}

/** 如果你希望将问题数据复制到子 SCIP (sub-SCIP) 中，则在此处复制用户数据 */
SCIP_RETCODE ProbDataTSP::scip_copy(
   SCIP*                 scip,               /**< SCIP 数据结构 */
   SCIP*                 sourcescip,         /**< 源 SCIP 主数据结构 */
   SCIP_HASHMAP*         varmap,             /**< 哈希表，存储源变量到对应目标变量的映射 */
   SCIP_HASHMAP*         consmap,            /**< 哈希表，存储源约束到对应目标约束的映射 */
   ObjProbData**         objprobdata,        /**< 用于存储复制后的问题数据对象的指针 */
   SCIP_Bool             global,             /**< 创建全局副本还是局部副本？ */
   SCIP_RESULT*          result              /**< 用于存储调用结果的指针 */
   )
{
   // 获取源问题数据及其图
   ProbDataTSP* sourceprobdatatsp;
   sourceprobdatatsp = dynamic_cast<ProbDataTSP*>(SCIPgetObjProbData(sourcescip));
   assert( sourceprobdatatsp != NULL );

   GRAPH* sourcegraph = sourceprobdatatsp->graph_;
   assert( sourcegraph != NULL );

   // 复制图
   GRAPH* graph = NULL;
   SCIP_CALL( copy_graph(&graph, sourcegraph) );

   // 复制并链接变量
   int m = graph->nedges;
   for(int e = 0; e < m; ++e)
   {
      SCIP_Bool success;
      GRAPHEDGE* edgeforw  = &(graph->edges[e]);
      GRAPHEDGE* edgebackw = &(graph->edges[e + m]);
      assert( sourcegraph->edges[e].var != NULL );

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcegraph->edges[e].var, &(edgeforw->var), varmap, consmap, global, &success) );
      SCIP_CALL( SCIPcaptureVar(scip, edgeforw->var) );
      assert(success);
      assert(edgeforw->var != NULL);

      // 反向平行的弧共享同一个变量
      edgebackw->var = edgeforw->var;
      SCIP_CALL( SCIPcaptureVar(scip, edgebackw->var) );
   }

   // 为目标问题数据分配内存
   ProbDataTSP* probdatatsp = new ProbDataTSP(graph);
   assert( probdatatsp != NULL );

   // 保存数据指针
   assert( objprobdata != NULL );
   *objprobdata = probdatatsp;

   // 图已被 ProbDataTSP(graph) 构造函数捕获（引用计数增加），因此这里释放局部指针
   release_graph(&graph);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** 用户问题数据的析构方法，用于释放原始用户数据（当原始问题被释放时调用） */
SCIP_RETCODE ProbDataTSP::scip_delorig(
   SCIP*                 scip                /**< SCIP 数据结构 */
   )
{
   for( int i = 0; i < graph_->nedges; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].back->var) );
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].var) );
   }
   release_graph(&graph_);

   return SCIP_OKAY;
}

/** 用户问题数据的析构方法，用于释放变换后的用户数据（当变换后的问题被释放时调用） */
SCIP_RETCODE ProbDataTSP::scip_deltrans(
   SCIP*                 scip                /**< SCIP 数据结构 */
   )
{
   for( int i = 0; i < graph_->nedges; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].back->var) );
      SCIP_CALL( SCIPreleaseVar(scip, &graph_->edges[i].var) );
   }
   release_graph(&graph_);

   return SCIP_OKAY;
}

/** 通过变换原始用户问题数据来创建变换后问题的用户数据
 *  （在问题被变换后调用）
 */
SCIP_RETCODE ProbDataTSP::scip_trans(
   SCIP*                 scip,               /**< SCIP 数据结构 */
   ObjProbData**         objprobdata,        /**< 用于存储变换后的问题数据对象的指针 */
   SCIP_Bool*            deleteobject        /**< 用于存储 SCIP 是否应在求解后删除该对象的标志 */
   )
{  /*lint --e{715}*/
   assert( objprobdata != NULL );
   assert( deleteobject != NULL );

   assert( graph_ != NULL );

   // 复制图
   GRAPH* transgraph = NULL;
   SCIP_CALL( copy_graph(&transgraph, graph_) );

   // 复制并链接变量（获取变换后的变量）
   int m = transgraph->nedges;
   for(int e = 0; e < m; ++e)
   {
      GRAPHEDGE* edgeforw  = &(transgraph->edges[e]);
      GRAPHEDGE* edgebackw = &(transgraph->edges[e + m]);
      assert( graph_->edges[e].var != NULL );

      SCIP_CALL( SCIPgetTransformedVar(scip, graph_->edges[e].var, &(edgeforw->var)) );
      SCIP_CALL( SCIPcaptureVar(scip, edgeforw->var) );

      edgebackw->var = edgeforw->var; // 反向平行的弧共享同一个变量
      assert( edgebackw->var != NULL );

      SCIP_CALL( SCIPcaptureVar(scip, edgebackw->var) );
   }

   // 为目标问题数据分配内存
   ProbDataTSP* transprobdatatsp = new ProbDataTSP(transgraph);
   assert( transprobdatatsp != NULL );

   // 保存数据指针
   assert( objprobdata != NULL );
   *objprobdata = transprobdatatsp;

   // 图已被 ProbDataTSP(graph) 构造函数捕获（引用计数增加），因此这里释放局部指针
   release_graph(&transgraph);

   *deleteobject = TRUE;

   return SCIP_OKAY;
}