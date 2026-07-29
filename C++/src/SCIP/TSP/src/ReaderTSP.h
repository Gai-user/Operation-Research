/**@file   ReaderTSP.h
 * @brief  TSP 数据文件的 C++ 文件读取器
 * @author Timo Berthold
 *
 * @details
 * 读取 TSPLIB 格式的 TSP 问题文件，解析节点坐标并构建完全图。
 *
 * ## TSPLIB 文件格式
 *
 * TSPLIB 是 TSP 问题的标准测试集格式。一个典型的文件结构如下：
 *
 * ~~~
 * NAME: gr17
 * TYPE: TSP
 * COMMENT: 17 cities
 * DIMENSION: 17
 * EDGE_WEIGHT_TYPE: EUC_2D
 * NODE_COORD_SECTION
 *   1  2.08  1.42
 *   2  1.68  1.90
 *   ...
 * EOF
 * ~~~
 *
 * ## 距离计算方式 (EDGE_WEIGHT_TYPE)
 *
 * | 类型 | 公式 | 说明 |
 * |------|------|------|
 * | EUC_2D | sqrt(dx² + dy²) | 二维欧几里得距离 |
 * | MAX_2D | max(\|dx\|, \|dy\|) | 切比雪夫距离 |
 * | MAN_2D | \|dx\| + \|dy\| | 曼哈顿距离 |
 * | GEO | 大圆距离 | 地理坐标（经纬度） |
 * | ATT | ceil(sqrt((dx²+dy²)/10)) | AT&T 伪欧几里得距离 |
 */

#pragma once

#include <string>
#include <vector>
#include <optional>
#include "objscip/objscip.h"
#include "scip/def.h"
#include "GomoryHuTree.h"
#include "scip/pub_fileio.h"

namespace tsp
{
    /**
     * @brief TSP 数据文件的 SCIP 文件读取器
     *
     * 负责从 TSPLIB 格式的文件中读取 TSP 实例，解析节点坐标，
     * 构建完全图（含半边模型），创建 SCIP 变量和初始约束。
     */
    class ReaderTSP : public scip::ObjReader
    {
    public:
        /**
         * @brief 构造读取器并注册参数
         *
         * 注册布尔参数 "reading/tspreader/round_lengths"：
         * 是否将边长四舍五入到最接近的整数。
         */
        explicit ReaderTSP(SCIP *scip)
            : scip::ObjReader(scip, "tspreader", "file reader for TSP files", "tsp")
        {
            (void)SCIPaddBoolParam(scip,
                                   "reading/tspreader/round_lengths",
                                   "should lengths of edges be rounded to nearest integer?",
                                   &round_lengths_, FALSE, TRUE, nullptr, nullptr);
        }

        virtual ~ReaderTSP() = default;

        /** @brief 读取器释放回调 */
        virtual SCIP_DECL_READERFREE(scip_free);
        /** @brief 文件读取回调：解析 TSPLIB 文件并构建问题 */
        virtual SCIP_DECL_READERREAD(scip_read);
        /** @brief 文件写入回调（暂未实现） */
        virtual SCIP_DECL_READERWRITE(scip_write);

    private:
        SCIP_Bool round_lengths_ = FALSE; ///< 是否对边长四舍五入

        /**
         * @brief 从文件中解析节点坐标
         *
         * 逐行读取 NODE_COORD_SECTION 中的节点编号和 (x, y) 坐标，
         * 填充到图的节点中。
         */
        SCIP_RETCODE getNodesFromFile(
            SCIP_FILE *file,
            std::vector<double> &x_coords,
            std::vector<double> &y_coords,
            Graph *graph);

        /**
         * @brief 验证 TSP 文件格式的有效性
         *
         * 检查：
         *   - 节点数量 > 0
         *   - TYPE = "TSP"
         *   - EDGE_WEIGHT_TYPE 为支持的类型之一
         *   - 图已成功初始化
         *
         * @return 错误信息的 optional，如果验证通过则为 nullopt
         */
        static std::optional<std::string> checkValid(
            const Graph *graph,
            const std::string &name,
            const std::string &type,
            const std::string &edge_weight_type,
            int num_nodes);
    };
}
