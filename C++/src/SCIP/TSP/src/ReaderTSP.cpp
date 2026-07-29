/**@file   ReaderTSP.cpp
 * @brief  TSP 数据文件的 C++ 文件读取器实现
 * @author Timo Berthold
 *
 * @details
 * 实现从 TSPLIB 格式文件读取 TSP 问题的完整流程：
 *
 * ## 读取流程
 *
 * 1. 打开文件，按行解析
 * 2. 解析文件头：NAME, TYPE, DIMENSION, EDGE_WEIGHT_TYPE
 * 3. 解析 NODE_COORD_SECTION：读取每个节点的 (x, y) 坐标
 * 4. 构建完全图：
 *    - 对每对节点 (i, j)，创建一对互为反向的半边
 *    - 计算边长（根据 EDGE_WEIGHT_TYPE）
 *    - 构建邻接表
 * 5. 创建 SCIP 模型：
 *    - 为每条无向边创建二元变量 x_e ∈ {0, 1}
 *    - 为每个节点创建度约束 Σ x_e = 2
 *    - 创建子环消除约束
 *
 * ## 半边索引约定
 *
 * 对于 n 个节点的完全图，有 m = n*(n-1)/2 条无向边，即 2m 条半边。
 * 半边在 edges_ vector 中的布局：
 *   - edges_[0..m-1]：所有"正向"半边（从小编号节点到大编号节点）
 *   - edges_[m..2m-1]：所有"反向"半边（从大编号节点到小编号节点）
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cctype>

#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "scip/pub_fileio.h"

#include "ReaderTSP.h"
#include "ProbDataTSP.h"
#include "ConshdlrSubtour.h"
#include "GomoryHuTree.h"

using namespace scip;

namespace tsp
{

/** 四舍五入到最近整数的宏（用于 ATT 距离的取整） */
#define NINT(x) (floor((x) + 0.5))

    /* ==========================================================================
     * Token 提取
     * ========================================================================== */

    /**
     * @brief 从 C 字符串中提取下一个空白分隔的 token
     *
     * 跳过开头的空白和冒号字符，收集连续的非空白非冒号字符作为 token。
     * 返回值后，str 指针指向下一段内容的起始位置。
     *
     * @param str 指向 C 字符串当前解析位置的指针（会被修改）
     * @return 提取的 token 字符串
     */
    static std::string getToken(const char *&str)
    {
        // 跳过空白字符和冒号（TSPLIB 中冒号是关键字分隔符）
        while (*str != '\0' && (std::isspace(static_cast<unsigned char>(*str)) || *str == ':'))
            ++str;
        // 收集连续的 token 字符
        std::string token;
        while (*str != '\0' && *str != ':' && !std::isspace(static_cast<unsigned char>(*str)))
        {
            token += *str;
            ++str;
        }
        // 跳过 token 后的空白和冒号
        while (*str != '\0' && (std::isspace(static_cast<unsigned char>(*str)) || *str == ':'))
            ++str;

        return token;
    }

    /* ==========================================================================
     * 节点解析
     * ========================================================================== */

    /**
     * @brief 从文件中逐行解析节点坐标
     *
     * 读取 NODE_COORD_SECTION 中的每一行，格式为：
     * ~~~
     *   <节点编号> <X坐标> <Y坐标>
     * ~~~
     *
     * 节点的编号必须与数组下标一致（从 1 开始编号，内部转为 0-based ID）。
     * 如果编号不匹配，发出警告并返回错误。
     *
     * @param file      已打开的 TSPLIB 文件
     * @param x_coords  输出：各节点的 X 坐标数组
     * @param y_coords  输出：各节点的 Y 坐标数组
     * @param graph     图（坐标写入 node.coord_x 和 node.coord_y）
     * @return SCIP_OKAY 成功，否则 SCIP_INVALIDDATA
     */
    SCIP_RETCODE ReaderTSP::getNodesFromFile(
        SCIP_FILE *file,
        std::vector<double> &x_coords,
        std::vector<double> &y_coords,
        Graph *graph)
    {
        char str[SCIP_MAXSTRLEN];
        int i = 0;
        int num_nodes = graph->numNodes();

        while (i < num_nodes && !SCIPfeof(file))
        {
            (void)SCIPfgets(str, SCIP_MAXSTRLEN, file);
            const char *s = str;

            // 解析节点编号
            int node_number;
            if (!SCIPstrToIntValue(str, &node_number, const_cast<char **>(&s)))
            {
                SCIPerrorMessage("Could not read node number:\n%s\n", str);
                return SCIP_INVALIDDATA;
            }

            // 解析 X 坐标
            if (!SCIPstrToRealValue(s, &x_coords[i], const_cast<char **>(&s)))
            {
                SCIPerrorMessage("Could not read x coordinate:\n%s\n", str);
                return SCIP_INVALIDDATA;
            }

            // 解析 Y 坐标
            if (!SCIPstrToRealValue(s, &y_coords[i], const_cast<char **>(&s)))
            {
                SCIPerrorMessage("Could not read y coordinate:\n%s\n", str);
                return SCIP_INVALIDDATA;
            }

            // 验证节点编号一致性
            auto &node = graph->node(i);
            node.id = i;
            if (node_number - 1 != i)
            {
                std::cout << "warning: nodenumber <" << node_number
                          << "> does not match its index in node list <" << i + 1
                          << ">. Node will get number " << i + 1
                          << " when naming variables and constraints!" << std::endl;
                return SCIP_INVALIDDATA;
            }

            // 写入坐标
            node.coord_x = x_coords[i];
            node.coord_y = y_coords[i];
            node.first_edge = nullptr;

            ++i;
        }

        assert(i == num_nodes);
        return SCIP_OKAY;
    }

    /* ==========================================================================
     * 文件格式验证
     * ========================================================================== */

    /**
     * @brief 验证 TSP 文件格式的有效性
     *
     * 进行以下检查：
     *   1. 节点数量 > 0
     *   2. TYPE 必须为 "TSP"
     *   3. EDGE_WEIGHT_TYPE 必须是支持的类型之一
     *   4. 图必须已成功初始化（NODE_COORD_SECTION 已解析）
     *
     * @return 验证失败时的错误信息；验证通过则返回 nullopt
     */
    std::optional<std::string> ReaderTSP::checkValid(
        const Graph *graph,
        const std::string &name,
        const std::string &type,
        const std::string &edge_weight_type,
        int num_nodes)
    {
        if (num_nodes < 1)
            return "parse error in file <" + name + "> dimension should be greater than 0";

        if (type != "TSP")
            return "parse error in file <" + name + "> type should be TSP";

        if (edge_weight_type != "EUC_2D" && edge_weight_type != "MAX_2D" &&
            edge_weight_type != "MAN_2D" && edge_weight_type != "GEO" &&
            edge_weight_type != "ATT")
            return "parse error in file <" + name + "> unknown weight type, "
                                                    "should be EUC_2D, MAX_2D, MAN_2D, ATT, or GEO";

        if (graph == nullptr)
            return "error while reading file <" + name + ">, graph is uninitialized. "
                                                         "Probably NODE_COORD_SECTION is missing";

        return std::nullopt;
    }

    /* ==========================================================================
     * 边长计算
     *
     * 支持 TSPLIB 中定义的五种距离计算方式。
     * ========================================================================== */

    /**
     * @brief 根据 EDGE_WEIGHT_TYPE 计算两点之间的距离
     *
     * ## 计算公式
     *
     * | 类型     | 公式                          | 说明                      |
     * |----------|-------------------------------|---------------------------|
     * | EUC_2D   | sqrt(dx² + dy²)               | 标准二维欧几里得距离        |
     * | MAX_2D   | max(|dx|, |dy|)              | 切比雪夫距离（棋盘距离）     |
     * | MAN_2D   | |dx| + |dy|                 | 曼哈顿距离（城市街区距离）   |
     * | ATT      | ceil(sqrt((dx²+dy²)/10))     | AT&T 伪欧几里得距离        |
     * | GEO      | 大圆距离公式                  | 地理坐标（经纬度）的球面距离 |
     *
     * @param edge_weight_type 距离计算方式标识符
     * @param x1, y1 第一个点的坐标
     * @param x2, y2 第二个点的坐标
     * @return 两点间的距离
     */
    static double computeEdgeLength(
        const std::string &edge_weight_type,
        double x1, double y1,
        double x2, double y2)
    {
        double dx = x1 - x2;
        double dy = y1 - y2;

        // 欧几里得距离：二维平面上两点间的直线距离
        if (edge_weight_type == "EUC_2D")
            return std::sqrt(dx * dx + dy * dy);

        // 切比雪夫距离：以最大坐标差为距离
        if (edge_weight_type == "MAX_2D")
            return std::max(std::abs(dx), std::abs(dy));

        // 曼哈顿距离：只能沿水平或垂直方向移动
        if (edge_weight_type == "MAN_2D")
            return std::abs(dx) + std::abs(dy);

        // AT&T 伪欧几里得距离：除以 10 后开方再向上取整
        if (edge_weight_type == "ATT")
            return std::ceil(std::sqrt((dx * dx + dy * dy) / 10.0));

        // GEO：地理坐标的大圆距离
        // 将 DDD.MM 格式的经纬度转换为弧度，然后使用球面余弦公式
        const double pi = 3.141592653589793;
        double coords[4] = {x1, y1, x2, y2};
        double rads[4];

        for (int k = 0; k < 4; ++k)
        {
            // 分离度数和分钟：整数部分为度，小数部分为分钟
            double deg = (coords[k] >= 0) ? std::floor(coords[k]) : std::ceil(coords[k]);
            double min = coords[k] - deg;
            // 转换为弧度：度 + 分钟/60 × π/180
            rads[k] = pi * (deg + 5.0 * min / 3.0) / 180.0;
        }

        // 球面余弦公式
        double q1 = std::cos(rads[1] - rads[3]); // cos(Δlat)
        double q2 = std::cos(rads[0] - rads[2]); // cos(Δlon)
        double q3 = std::cos(rads[0] + rads[2]); // cos(sum lon)

        // 地球半径 6378.388 km × arccos(距离余弦)，结果向下取整 + 1
        return std::floor(6378.388 * std::acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    }

    /* ==========================================================================
     * 构建完全图
     * ========================================================================== */

    /**
     * @brief 根据节点坐标和距离方式构建带权完全图
     *
     * 对每对节点 (i, j), i < j，创建一对互为反向的半边，计算边长，
     * 并链接到各自的邻接表。
     *
     * 半边索引：使用 edge_idx 计数器从 0 递增。
     * 正向半边存储在 edges[edge_idx]，反向半边存储在 edges[edge_idx + num_edges]。
     *
     * @param graph             已分配节点但未连接边的图
     * @param x_coords          节点的 X 坐标数组
     * @param y_coords          节点的 Y 坐标数组
     * @param edge_weight_type  距离计算方式
     * @param round_lengths     是否对边长进行四舍五入
     */
    static void buildDistanceGraph(
        Graph &graph,
        const std::vector<double> &x_coords,
        const std::vector<double> &y_coords,
        const std::string &edge_weight_type,
        bool round_lengths)
    {
        int num_nodes = graph.numNodes();
        int num_edges = graph.numEdges();

        int edge_idx = 0;
        for (int i = 0; i < num_nodes; ++i)
        {
            auto &start_node = graph.node(i);
            for (int j = i + 1; j < num_nodes; ++j)
            {
                auto &end_node = graph.node(j);

                // 获取半边
                auto &forward_edge = graph.edge(edge_idx);
                auto &reverse_edge = graph.edge(edge_idx + num_edges);

                // 设置目标节点
                forward_edge.target = &end_node;
                reverse_edge.target = &start_node;

                forward_edge.setReverseEdge(&reverse_edge);

                // 计算边长
                double length = computeEdgeLength(edge_weight_type,
                                                  x_coords[start_node.id], y_coords[start_node.id],
                                                  x_coords[end_node.id], y_coords[end_node.id]);

                if (round_lengths)
                {
                    length = NINT(length);
                }

                // 设置边长
                forward_edge.length = length;
                reverse_edge.length = length;

                // 链接到邻接表头部
                forward_edge.linkTo(&start_node);
                reverse_edge.linkTo(&end_node);

                ++edge_idx;
            }
        }
    }

    /**
     * @brief 为图的每条无向边创建 SCIP 二元决策变量
     *
     * 每条无向边创建一个二元变量 x_e ∈ {0, 1}：
     *   - 下界 = 0，上界 = 1
     *   - 目标函数系数 = 边长 (length)
     *   - 变量类型 = BINARY
     *   - 变量名格式：x_e_<node1>-<node2>
     *
     * 两条半边（正向和反向）通过 setVar() 共享同一个变量引用。
     *
     * @param scip   SCIP 实例
     * @param graph  图
     * @return SCIP_OKAY 成功，否则错误码
     */
    static SCIP_RETCODE createSCIPVariables(SCIP *scip, Graph *graph)
    {
        int num_edges = graph->numEdges();

        for (int i = 0; i < num_edges; ++i)
        {
            auto &edge = graph->edge(i);
            SCIP_VAR *var;

            // 生成变量名：x_e_<start>-<end>
            std::ostringstream varname;
            varname << "x_e_" << edge.reverse->target->id + 1
                    << "-" << edge.target->id + 1;

            // 创建二元变量：
            //   下界 0.0，上界 1.0，目标系数 = 边长
            //   BINARY 类型，需要在解和 LP 中都出现
            SCIP_CALL(SCIPcreateVar(scip, &var, varname.str().c_str(),
                                    0.0, 1.0, edge.length,
                                    SCIP_VARTYPE_BINARY, TRUE, FALSE,
                                    nullptr, nullptr, nullptr, nullptr, nullptr));

            // 将变量添加到SCIP实例
            SCIP_CALL(SCIPaddVar(scip, var));

            // 设置边的变量
            SCIP_CALL(edge.setVar(scip, var));

            SCIP_CALL(SCIPreleaseVar(scip, &var));
        }

        return SCIP_OKAY;
    }

    /**
     * @brief 为每个节点创建度约束
     *
     * 在 TSP 中，每个节点恰好有两条邻接边被选中（一条进入，一条离开）：
     * ~~~
     *   Σ_{e ∈ δ(v)} x_e = 2   对于所有节点 v
     * ~~~
     *
     * 其中 δ(v) 是与节点 v 关联的所有边的集合。
     *
     * 度约束是 TSP 的基本约束。它保证了每个节点的连通性，
     * 但不足以防止子环出现（多个不相交的环也满足度约束）。
     * 子环消除需要额外的约束处理器。
     *
     * @param scip   SCIP 实例
     * @param graph  图
     * @return SCIP_OKAY 成功，否则错误码
     */
    static SCIP_RETCODE createInitialDegreeConstraints(SCIP *scip, Graph *graph)
    {
        int num_nodes = graph->numNodes();

        if (num_nodes < 2)
        {
            return SCIP_OKAY;
        }

        for (int i = 0; i < num_nodes; ++i)
        {
            auto &node = graph->node(i);
            SCIP_CONS *cons;
            std::ostringstream consname;
            consname << "deg_con_v" << node.id + 1;

            // 创建线性约束：左侧 = 2.0, 右侧 = 2.0（等式约束）
            // 参数：左侧 = 2.0, 右侧 = 2.0 → 等式
            SCIP_CALL(SCIPcreateConsLinear(scip, &cons, consname.str().c_str(),
                                           0, nullptr, nullptr, 2.0, 2.0,
                                           TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

            // 累加该节点所有邻接边关联的变量
            for (auto *edge = node.first_edge; edge != nullptr; edge = edge->next)
            {
                SCIP_CALL(SCIPaddCoefLinear(scip, cons, edge->scip_var, 1.0));
            }

            // 将约束添加到SCIP
            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));
        }

        return SCIP_OKAY;
    }

    /* ==========================================================================
     * SCIP reader 回调函数
     * ========================================================================== */

    /**
     * @brief 读取器释放回调（无动态分配资源需要释放）
     */
    SCIP_DECL_READERFREE(ReaderTSP::scip_free)
    {
        return SCIP_OKAY;
    }

    /**
     * @brief 文件读取回调：完整解析 TSPLIB 文件并构建 SCIP 模型
     *
     * ## 执行流程
     *
     * 1. 打开文件
     * 2. 解析文件头（NAME, TYPE, DIMENSION, EDGE_WEIGHT_TYPE）
     * 3. 解析 NODE_COORD_SECTION，创建图
     * 4. 验证文件格式
     * 5. 构建完全图（连接所有节点对，计算边长）
     * 6. 创建问题数据对象（ProbDataTSP）
     * 7. 创建 SCIP 二元变量（每条边一个）
     * 8. 创建度约束（每个节点 Σx_e = 2）
     * 9. 创建子环消除约束
     *
     * @param scip     SCIP 实例
     * @param reader   读取器对象
     * @param filename 要读取的文件路径
     * @param result   输出：读取结果
     */
    SCIP_DECL_READERREAD(ReaderTSP::scip_read)
    {
        SCIP_RETCODE retcode = SCIP_OKAY;
        *result = SCIP_DIDNOTRUN;

        //  打开文件
        SCIP_FILE *file = SCIPfopen(filename, "r");
        if (!file)
        {
            return SCIP_READERROR;
        }

        // 初始化默认值
        std::string name = "MY_OWN_LITTLE_TSP";
        std::string type = "TSP";
        std::string edge_weight_type = "EUC_2D";
        int num_nodes = 0;
        int num_edges = 0;
        std::unique_ptr<Graph> graph;
        std::vector<double> x_coords;
        std::vector<double> y_coords;

        // 逐行解析文件
        char str[SCIP_MAXSTRLEN];
        (void)SCIPfgets(str, SCIP_MAXSTRLEN, file);
        const char *s = str;
        std::string token = getToken(s);

        while (!SCIPfeof(file))
        {
            if (token == "NAME")
            {
                name = getToken(s);
            }
            else if (token == "TYPE")
            {
                type = getToken(s);
            }
            else if (token == "DIMENSION")
            {
                if (!SCIPstrToIntValue(s, &num_nodes, const_cast<char **>(&s)))
                {
                    SCIPerrorMessage("Could not read number of nodes:\n%s\n", s);
                    return SCIP_INVALIDDATA;
                }
                // 完全图的半边数量 = n * (n - 1)
                num_edges = num_nodes * (num_nodes - 1);
            }
            else if (token == "EDGE_WEIGHT_TYPE")
            {
                edge_weight_type = getToken(s);
            }
            else if (token == "NODE_COORD_SECTION" || token == "DISPLAY_DATA_SECTION")
            {
                // 解析节点坐标
                if (num_nodes < 1)
                {
                    retcode = SCIP_READERROR;
                    break;
                }

                // 创建图
                graph = Graph::create(num_nodes, num_edges);

                if (!graph)
                {
                    retcode = SCIP_NOMEMORY;
                    break;
                }

                x_coords.resize(num_nodes);
                y_coords.resize(num_nodes);
                SCIP_CALL(getNodesFromFile(file, x_coords, y_coords, graph.get()));
            }
            else if (token == "COMMENT:" || token == "COMMENT" ||
                     token == "DISPLAY_DATA_TYPE" || token == "DISPLAY_DATA_TYPE:")
            {
                // 跳过注释行和显示数据类型行
            }
            else if (token == "EOF")
            {
                break;
            }
            else if (token == "")
            {
                // 跳过空行
            }
            else
            {
                std::cout << "parse error in file <" << name << "> unknown keyword <"
                          << token << ">" << std::endl;
                return SCIP_READERROR;
            }

            (void)SCIPfgets(str, SCIP_MAXSTRLEN, file);
            s = str;
            token = getToken(s);
        }

        SCIPfclose(file);

        // 验证文件格式
        auto err = checkValid(graph.get(), name, type, edge_weight_type, num_nodes);
        if (err)
        {
            std::cout << *err << std::endl;
            retcode = SCIP_READERROR;
        }

        assert(graph != nullptr);

        if (retcode != SCIP_OKAY)
        {
            return retcode;
        }

        // 构建完全图
        buildDistanceGraph(*graph, x_coords, y_coords, edge_weight_type, round_lengths_);

        // 创建问题数据对象
        SCIP_CALL(SCIPcreateObjProb(scip, name.c_str(), new ProbDataTSP(std::move(graph)), true));

        // 获取图的引用
        ProbDataTSP *probdata = dynamic_cast<ProbDataTSP *>(SCIPgetObjProbData(scip));
        assert(probdata != nullptr);
        Graph *g = probdata->getGraph();

        // 创建变量、约束
        SCIP_CALL(createSCIPVariables(scip, g));
        SCIP_CALL(createInitialDegreeConstraints(scip, g));

        // 创建子环消除约束

        SCIP_CONS *cons;
        SCIP_CALL(SCIPcreateConsSubtour(scip, &cons, "subtour", g,
                                        FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));
        SCIP_CALL(SCIPaddCons(scip, cons));
        SCIP_CALL(SCIPreleaseCons(scip, &cons));

        *result = SCIP_SUCCESS;
        return SCIP_OKAY;
    }

    /**
     * @brief 文件写入回调（暂未实现）
     */
    SCIP_DECL_READERWRITE(ReaderTSP::scip_write)
    {
        *result = SCIP_DIDNOTRUN;
        return SCIP_OKAY;
    }
} // end tsp