/**@file
 * @brief  VRP 定价器示例的主文件
 * @author Andreas Bley
 * @author Marc Pfetsch
 *
 * 我们希望在图 \f$G = (V,E)\f$ 上求解车辆路径问题（VRP），其中
 * \f$V = J \cup \{d\}\f$，d 表示仓库（depot），距离由
 * 长度函数 \f$l_e: E \rightarrow R_{\geq 0}\f$ 给出。
 *
 * 考虑如下的混合整数规划（MIP）模型：
 *
 * \f[
 *  \begin{array}[t]{rll}
 *    \min &  \displaystyle \sum_{e \in E} l_e y_e \\
 *         & & \\
 *   s.t.  & -y_e + \sum_{t \in T_k} a^t_e x_t  \leq 0, &  \forall e \in E\\
 *         &  \displaystyle \sum_{t \in T_k} a^t_j x_t = 1, &  \forall j \in J \\
 *         &  y(\delta(j)) = 2, &  \forall j \in J \\
 *         &  y_e \in \{0,1,2\},  & \forall e \in E \\
 *         &  x_t  \in [0,1], & \forall t \in T_k
 *  \end{array}
 * \f]
 *
 * 其中，\f$T_k\f$ 是访问最多 k 个客户的路线集合，
 * 允许重复访问客户，\f$a^t_e\f$（\f$a^t_j\f$）表示在 \f$t \in T_k\f$ 中
 * 边 e（节点 j）被经过的次数。
 *
 * 相关示例和文件格式详见：https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-instances/。
 */

/* 标准库 */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

/* SCIP 相关 */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "PricerVRP.h"

using namespace std;
using namespace scip;

/* 读取问题*/
static int readProblem(
    const char *filename,     /**< filename */
    int &numNodes,           /**< number of nodes in instance */
    int &capacity,            /**< capacity in instance */
    vector<int> &demand,      /**< array of demands of instance */
    vector<vector<int>> &dist /**< distances between nodes */
)
{
    const string DIMENSION = "DIMENSION";
    const string DEMAND_SECTION = "DEMAND_SECTION";
    const string DEPOT_SECTION = "DEPOT_SECTION";
    const string EDGE_WEIGHT_TYPE = "EDGE_WEIGHT_TYPE";
    const string EUC_2D = "EUC_2D";
    const string EXPLICIT = "EXPLICIT";
    const string LOWER_DIAG_ROW = "LOWER_DIAG_ROW";
    const string EDGE_WEIGHT_FORMAT = "EDGE_WEIGHT_FORMAT";
    const string EDGE_WEIGHT_SECTION = "EDGE_WEIGHT_SECTION";
    const string NODE_COORD_SECTION = "NODE_COORD_SECTION";
    const string CAPACITY = "CAPACITY";

    ifstream file(filename);

    numNodes = -1;
    capacity = -1;

    if (!file)
    {
        cerr << "Cannot open file" << filename << endl;
        return 1;
    }

    string edge_weight_type = "";
    string edge_weight_format = "";
    vector<int> x;
    vector<int> y;

    while (file)
    {
        string key;
        string dummy;
        file >> key;

        if (key == DIMENSION)
        {
            file >> dummy;
            file >> numNodes;

            assert(numNodes >= 0);
            demand.resize(numNodes, 0);
            dist.resize(numNodes);
            for (int i = 0; i < numNodes; ++i)
            {
                dist[i].resize(i, 0);
            }
        }

        if (key == CAPACITY)
        {
            file >> dummy;
            file >> capacity;
        }
        else if (key == EDGE_WEIGHT_TYPE)
        {
            file >> dummy;
            file >> edge_weight_type;

            if (edge_weight_type != EUC_2D && edge_weight_type != EXPLICIT)
            {
                cerr << "Wrong " << EDGE_WEIGHT_TYPE << " " << edge_weight_type << endl;
                return 1;
            }
            if (edge_weight_type == EUC_2D)
            {
                assert(numNodes >= 0);
                x.resize(numNodes, 0);
                y.resize(numNodes, 0);
            }
        }
        else if (key == EDGE_WEIGHT_FORMAT)
        {
            file >> dummy;
            file >> edge_weight_format;
        }
        else if (key == EDGE_WEIGHT_FORMAT + ":")
        {
            file >> edge_weight_format;
        }
        else if (key == EDGE_WEIGHT_SECTION)
        {
            if (edge_weight_type != EXPLICIT || edge_weight_format != LOWER_DIAG_ROW)
            {
                cerr << "Error. Unsupported edge length type." << endl;
                return 1;
            }
            assert(numNodes >= 0);
            for (int i = 0; i < numNodes; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    int l;
                    file >> l;
                    dist[i][j] = l;
                }
            }
        }
        else if (key == NODE_COORD_SECTION)
        {
            if (edge_weight_type != EUC_2D)
            {
                cerr << "Error. Data file contains " << EDGE_WEIGHT_TYPE << " " << edge_weight_type << " and " << NODE_COORD_SECTION << endl;
                return 1;
            }
            assert(numNodes >= 0);
            for (int i = 0; i < numNodes; ++i)
            {
                int j, xi, yi;
                file >> j;
                file >> xi;
                file >> yi;
                if (j != i + 1)
                {
                    cerr << "Error reading " << NODE_COORD_SECTION << endl;
                    return 1;
                }
                x[i] = xi; /*lint !e732 !e747*/
                y[i] = yi; /*lint !e732 !e747*/
            }
            for (int i = 0; i < numNodes; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    int dx = x[i] - x[j];                                    /*lint !e732 !e747 !e864*/
                    int dy = y[i] - y[j];                                    /*lint !e732 !e747 !e864*/
                    dist[i][j] = int(sqrt((double)dx * dx + dy * dy) + 0.5); /*lint !e732 !e747 !e790*/
                }
            }
        }
        else if (key == DEMAND_SECTION)
        {
            assert(numNodes >= 0);
            for (int i = 0; i < numNodes; ++i)
            {
                int j, d;
                file >> j;
                file >> d;
                if (j != i + 1)
                {
                    cerr << "Error reading " << DEMAND_SECTION << endl;
                    return 1;
                }
                demand[i] = d; /*lint !e732 !e747*/
            }
        }
        else if (key == DEPOT_SECTION)
        {
            for (int i = 0; i != -1;)
            {
                file >> i;
                if (i != -1 && i != 1)
                {
                    cerr << "Error: This file specifies other depots than 1." << endl;
                    return 1;
                }
            }
        }
        else
        {
            (void)getline(file, dummy);
        }
    }
    return 0;
}

static SCIP_RETCODE runSCIP(int argc, char **argv)
{
    SCIP *scip = nullptr;

    cout << "Solving the vehicle routing problem using SCIP" << endl;
    cout << "Implemented bu GaiGuoNing" << endl;

    if (argc != 2 && argc != 3)
    {
        cerr << "Usage: vrp [-h] datafile" << endl;
        cerr << "Options:" << endl;
        cerr << "-h Uses hop limit instead of capacity limit for tours" << endl;
        return SCIP_INVALIDDATA;
    }

    /* 问题数据 */

    const char *VRP_PRICER_NAME = "VRP_Pricer";

    vector<vector<int>> dist;
    vector<int> demand;
    int capacity;
    int numNodes;

    // TODO: implement readProblem
    if (readProblem(argv[argc - 1], numNodes, capacity, demand, dist))
    {
        cerr << "Error reading data file " << argv[argc - 1] << endl;
        return SCIP_READERROR;
    }

    assert(numNodes >= 0);
    assert(capacity >= 0);

    cout << "Number of nodes: " << numNodes << endl;

    if (argc == 3)
    {
        if (string("-h") != argv[1])
        {
            cerr << "Unkonw option " << argv[2] << endl;
            return SCIP_PARAMETERUNKNOWN;
        }

        int totalDemand = 0;
        for (int i = 1; i < numNodes; ++i)
        {
            totalDemand += demand[i];
        }

        if (totalDemand == 0.0)
        {
            cerr << "Total demand is zero!" << endl;
            return SCIP_INVALIDDATA;
        }

        capacity = (numNodes - 1) * capacity / totalDemand;
        // clear the vector and initialize all elements to 1
        demand.assign(numNodes, 1);
        demand[0] = 0;
        cout << "Max customers per tour: " << capacity << endl;
    }
    else
    {
        cout << "Max demand per tour: " << capacity << endl;
    }

    /* 设置 SCIP */

    /* 初始化 SCIP 环境*/
    SCIP_CALL(SCIPcreate(&scip));

    /* 版本信息 */
    SCIPprintVersion(scip, nullptr);
    SCIPinfoMessage(scip, nullptr, "\n");

    /* 注册默认插件 */
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));

    /*设置求解日志打印程度 */
    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 5));

    /* 创建一个空的问题 */
    SCIP_CALL(SCIPcreateProb(scip, "VRP", 0, 0, 0, 0, 0, 0, 0));

    /* add arc-routing variables */
    char var_name[255];
    vector<vector<SCIP_VAR *>> arcVar(numNodes);

    for (int i = 0; i < numNodes; ++i)
    {
        arcVar[i].resize(i, (SCIP_VAR *)nullptr);
        for (int j = 0; j < i; ++j)
        {
            SCIP_VAR *var;
            (void)SCIPsnprintf(var_name, 255, "E%d_%d", i, j);

            SCIP_CALL(SCIPcreateVar(scip,
                                    &var,
                                    var_name,
                                    0.0,
                                    2.2,
                                    dist[i][j],
                                    SCIP_VARTYPE_INTEGER,
                                    TRUE,
                                    FALSE,
                                    NULL, NULL, NULL, NULL, NULL));
            SCIP_CALL(SCIPaddVar(scip, var));
            arcVar[i][j] = var;
        }
    }

    /* add arc-routing - tour constraints*/
    char con_name[255];
    vector<vector<SCIP_CONS *>> arcCon(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
        arcCon[i].resize(i, (SCIP_CONS *)nullptr);
        for (int j = 0; j < i; ++j)
        {
            SCIP_CONS *con;
            (void)SCIPsnprintf(con_name, 255, "A%d_%d", i, j);
            SCIP_VAR *idx = arcVar[i][j];
            SCIP_Real coeff = -1;
            SCIP_CALL(SCIPcreateConsLinear(scip, &con, con_name, 1, &idx, &coeff,
                                           -SCIPinfinity(scip), /* lhs */
                                           0.0,                 /* rhs */
                                           TRUE,                /* initial */
                                           FALSE,               /* separate */
                                           TRUE,                /* enforce */
                                           TRUE,                /* check */
                                           TRUE,                /* propagate */
                                           FALSE,               /* local */
                                           TRUE,                /* modifiable */
                                           FALSE,               /* dynamic */
                                           FALSE,               /* removable */
                                           FALSE));             /* stickingatnode */
            SCIP_CALL(SCIPaddCons(scip, con));
            arcCon[i][j] = con;
        }
    }

    /* add arc-routing - degree constraints */
    for (int i = 1; i < numNodes; ++i)
    {
        SCIP_CONS *con;
        (void)SCIPsnprintf(con_name, 255, "D%d", i);
        SCIP_CALL(SCIPcreateConsLinear(scip, &con, con_name, 0, 0, 0,
                                       2.0,     /* lhs */
                                       2.0,     /* rhs */
                                       TRUE,    /* initial */
                                       FALSE,   /* separate */
                                       TRUE,    /* enforce */
                                       TRUE,    /* check */
                                       TRUE,    /* propagate */
                                       FALSE,   /* local */
                                       FALSE,   /* modifiable */
                                       FALSE,   /* dynamic */
                                       FALSE,   /* removable */
                                       FALSE)); /* stickingatnode */
        SCIP_CALL(SCIPaddCons(scip, con));
        for (int j = 0; j < numNodes; ++j)
        {
            if (j == i)
            {
                continue;
            }

            SCIP_CALL(SCIPaddCoefLinear(scip, con, i > j ? arcVar[i][j] : arcVar[j][i], 1.0));
        }
        SCIP_CALL(SCIPreleaseCons(scip, &con));
    }

    /* add set packing constraints (Node 0 is the depot) */
    vector<SCIP_CONS *> partCon(numNodes, (SCIP_CONS *)nullptr);
    for (int i = 1; i < numNodes; ++i)
    {
        SCIP_CONS *con = nullptr;
        (void)SCIPsnprintf(con_name, 255, "C%d", i);
        SCIP_CALL(SCIPcreateConsLinear(scip, &con, con_name, 0, NULL, NULL,
                                       1.0,                /* lhs */
                                       SCIPinfinity(scip), /* rhs */
                                       TRUE,               /* initial */
                                       FALSE,              /* separate */
                                       TRUE,               /* enforce */
                                       TRUE,               /* check */
                                       TRUE,               /* propagate */
                                       FALSE,              /* local */
                                       TRUE,               /* modifiable */
                                       FALSE,              /* dynamic */
                                       FALSE,              /* removable */
                                       FALSE /* stickingatnode */));
        SCIP_CALL(SCIPaddCons(scip, con));
        partCon[i] = con; /*lint !e732 !e747*/
    }

    /* 注册自定义的VRP定价类*/
    ObjPricerVRP *vrpPricerPtr = new ObjPricerVRP(scip, VRP_PRICER_NAME, numNodes, capacity, demand, dist,
                                                  arcVar, arcCon, partCon);

    SCIP_CALL(SCIPincludeObjPricer(scip, vrpPricerPtr, true));

    SCIP_CALL(SCIPactivatePricer(scip, SCIPfindPricer(scip, VRP_PRICER_NAME)));

    /* solve problem */
    SCIP_CALL(SCIPsolve(scip));

    /* 打印求解信息 */
    SCIP_CALL(SCIPprintStatistics(scip, NULL));
    SCIP_CALL(SCIPprintBestSol(scip, NULL, FALSE));

    /* 释放 变量以及约束 */
    for (int i = 0; i < numNodes; ++i)
    {
        if (i > 0)
        {
            SCIP_CALL(SCIPreleaseCons(scip, &partCon[i]));
        }
        for (int j = 0; j < i; ++j)
        {
            SCIP_CALL(SCIPreleaseCons(scip, &arcCon[i][j]));
            SCIP_CALL(SCIPreleaseVar(scip, &arcVar[i][j]));
        }
    }
    SCIP_CALL(SCIPfree(&scip));

    BMScheckEmptyMemory();

    return SCIP_OKAY;
}

int main(int argc, char **argv)
{
    return runSCIP(argc, argv) != SCIP_OKAY ? 1 : 0;
}