/**@file   GomoryHuTree.cpp
 * @brief  无向图中全局割（global cuts）的生成器
 * @author Georg Skorobohatyj
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdio.h>
#include <assert.h>

#include "objscip/objscip.h"
#include "GomoryHuTree.h"

/** 用于数值比较的 epsilon 值 */
#define  EPS  1.0E-10

/* 静态变量 */
static GRAPHNODE** active;
static long* number;
static long max_dist;
static long bound;
static SCIP_Bool co_check;


/** 创建一个图 */
SCIP_Bool create_graph(
   int                   n,                  /**< 节点数量 */
   int                   m,                  /**< 边数量（半边缘/halfedges 的总数） */
   GRAPH**               gr                  /**< 用于存储图指针的指针 */
   )
{
   assert( gr != NULL );

   BMSallocMemory(gr);
   if( *gr == NULL )
      return FALSE;

   BMSallocMemoryArray(&(*gr)->nodes, n);
   if( (*gr)->nodes == NULL )
   {
      BMSfreeMemory(gr);
      return FALSE;
   }

   BMSallocMemoryArray(&(*gr)->edges, m);
   if( (*gr)->edges == NULL )
   {
      BMSfreeMemoryArray(&(*gr)->nodes);
      BMSfreeMemory(gr);
      return FALSE;
   }
   (*gr)->nuses = 1;
   (*gr)->nnodes = n;
   (*gr)->nedges = m/2;
   (*gr)->nedgesnonzero = m/2;

   return TRUE;
}

/** 释放一个图 */
static
void free_graph(
   GRAPH**               gr                  /**< 指向图的指针 */
   )
{
   assert(gr != NULL);
   assert(*gr != NULL);
   assert((*gr)->nuses == 0);

   BMSfreeMemory(&(*gr)->nodes);
   BMSfreeMemory(&(*gr)->edges);
   BMSfreeMemory(gr);
}

/** 引用计数增加（捕获图） */
void capture_graph(
   GRAPH*                gr                  /**< 图 */
   )
{
   assert(gr != NULL);

   ++gr->nuses;
}

/** 引用计数减少并释放图（如果计数归零） */
void release_graph(
   GRAPH**               gr                  /**< 指向图的指针 */
   )
{
   assert(gr != NULL);
   assert(*gr != NULL);

   --(*gr)->nuses;

   if( (*gr)->nuses == 0 )
      free_graph(gr);
   *gr = NULL;
}
