/**@file   ProbDataTSP.h
 * @brief  TSP 问题的 C++ 问题数据类
 * @author GaiGuoNing
 */

#ifndef __TSPPROBDATA_H__
#define __TSPPROBDATA_H__

#include "objscip/objscip.h"
#include "GomoryHuTree.h"

namespace tsp
{

/** TSP 问题的 SCIP 用户自定义问题数据 */
class ProbDataTSP : public scip::ObjProbData /*lint --e{3713}*/
{
   GRAPH*                graph_;             /**< 图数据 */

public:

   /** 默认构造函数 */
   ProbDataTSP(
      GRAPH*             g                   /**< 图数据 */
      )
      : graph_(g)
   {
      capture_graph(graph_);
   }

   /** 析构函数 */
   virtual ~ProbDataTSP()
   {
      if( graph_ != NULL )
         release_graph(&graph_); /*lint !e1551*/
   }

   /** 如果你希望将问题数据复制到子 SCIP (sub-SCIP) 中，则在此处复制用户数据 */
   virtual SCIP_RETCODE scip_copy(
      SCIP*              scip,               /**< SCIP 数据结构 */
      SCIP*              sourcescip,         /**< 源 SCIP 主数据结构 */
      SCIP_HASHMAP*      varmap,             /**< 哈希表，存储源变量到对应目标变量的映射 */
      SCIP_HASHMAP*      consmap,            /**< 哈希表，存储源约束到对应目标约束的映射 */
      ObjProbData**      objprobdata,        /**< 用于存储复制后的问题数据对象的指针 */
      SCIP_Bool          global,             /**< 创建全局副本还是局部副本？ */
      SCIP_RESULT*       result              /**< 用于存储调用结果的指针 */
      );

   /** 用户问题数据的析构方法，用于释放原始用户数据（当原始问题被释放时调用） */
   virtual SCIP_RETCODE scip_delorig(
      SCIP*              scip                /**< SCIP 数据结构 */
      );

   /** 用户问题数据的析构方法，用于释放变换后的用户数据（当变换后的问题被释放时调用） */
   virtual SCIP_RETCODE scip_deltrans(
      SCIP*              scip                /**< SCIP 数据结构 */
      );

   /** 通过变换原始用户问题数据来创建变换后问题的用户数据
    *  （在问题被变换后调用）
    */
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP 数据结构 */
      ObjProbData**      objprobdata,        /**< 用于存储变换后的问题数据对象的指针 */
      SCIP_Bool*         deleteobject        /**< 用于存储 SCIP 是否应在求解后删除该对象的标志 */
      );

   /* 获取图 */
   GRAPH* getGraph()
   {
      return graph_; /*lint !e1535*//*lint !e1806*/
   }

};/*lint !e1712*/

} /* namespace tsp */

#endif