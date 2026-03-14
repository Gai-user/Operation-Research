/**@file   ReaderTSP.h
 * @brief  用于 TSP 数据文件的 C++ 文件读取器
 * @author GaiGuoNing
 */

#ifndef __TSPREADER_H__
#define __TSPREADER_H__


#include <string>
#include "objscip/objscip.h"
#include "scip/def.h"
#include "GomoryHuTree.h"   // 用于 GRAPH 结构体
#include "scip/pub_fileio.h"


namespace tsp
{

/** 用于 TSP 数据文件的 SCIP 文件读取器 */
class ReaderTSP : public scip::ObjReader /*lint --e{3713}*/
{
public:

   SCIP_Bool round_lengths_;

   /** 默认构造函数 */
   ReaderTSP(SCIP* scip)
      : scip::ObjReader(scip, "tspreader", "file reader for TSP files", "tsp")
   {
      /* 添加 TSP 读取器参数 */
      (void) SCIPaddBoolParam(scip,
         "reading/tspreader/round_lengths", "是否应将边的长度四舍五入到最近的整数？",
         &round_lengths_, FALSE, TRUE, NULL, NULL);
   }

   /** 析构函数 */
   virtual ~ReaderTSP()
   {
   }

   /** 文件读取器的析构方法，用于释放用户数据（当 SCIP 退出时调用） */
   virtual SCIP_DECL_READERFREE(scip_free);

   /** 读取器的问题读取方法
    *
    *  *result 的可能返回值：
    *  - SCIP_SUCCESS    : 读取器正确读取了文件并创建了相应的问题
    *  - SCIP_DIDNOTRUN  : 该读取器不负责处理给定的输入文件
    *
    *  如果读取器在输入文件中检测到错误，它应返回 RETCODE SCIP_READERR 或 SCIP_NOFILE。
    */
   virtual SCIP_DECL_READERREAD(scip_read);

   /** 读取器的问题写入方法；注意：如果参数 "genericnames" 为 TRUE，则
    *  SCIP 已经将所有变量和约束的名称设置为通用名称；因此，此
    *  方法应始终使用 SCIPvarGetName() 和 SCIPconsGetName()；
    *
    *  *result 的可能返回值：
    *  - SCIP_SUCCESS    : 读取器正确读取了文件并创建了相应的问题
    *  - SCIP_DIDNOTRUN  : 该读取器不负责处理给定的输入文件
    *
    *  如果读取器在写入文件流时检测到错误，它应返回
    *  RETCODE SCIP_WRITEERROR。
    */
   virtual SCIP_DECL_READERWRITE(scip_write);

private:

   /** 解析节点列表 */
   SCIP_RETCODE getNodesFromFile(
      SCIP_FILE*         file,               /**< 包含要提取数据的文件 */
      double*            x_coords,           /**< 用于填充节点 x 坐标的 double 数组 */
      double*            y_coords,           /**< 用于填充节点 y 坐标的 double 数组 */
      GRAPH*             graph               /**< 将由节点生成的图 */
      );

   /** 断言文件格式正确且所有内容均已正确设置的方法 */
   bool checkValid(
      GRAPH*             graph,              /**< 构建的图，不应为 NULL */
      const std::string& name,               /**< 文件名 */
      const std::string& type,               /**< 问题类型，应为 "TSP" */
      const std::string& edgeweighttype,     /**< 边权重类型，应为 "EUC_2D", "MAX_2D", "MAN_2D",
                                              *   "ATT", 或 "GEO" */
      int                nnodes              /**< 问题的维度（节点数），至少应为 1 */
      );

   /** 将一个变量添加到两条半边缘（halfedges），并捕获它以供在图中使用 */
   SCIP_RETCODE addVarToEdges(
      SCIP*              scip,               /**< SCIP 数据结构 */
      GRAPHEDGE*         edge,               /**< 图的一条边 */
      SCIP_VAR*          var                 /**< 对应该边的变量 */
      );

};/*lint !e1712*/

} /* namespace tsp */

#endif