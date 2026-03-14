#include <iostream>
#include <string>
#include <sstream>

#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "scip/pub_fileio.h"
#include <math.h>

#include "ReaderTSP.h"
#include "ProbData.h"
#include "GomoryHuTree.h"

using namespace tsp;
using namespace scip;
using namespace std;

#define NINT(x) (floor(x + 0.5))

/** get token */
string getToken(char *&str)
{
    string token;

    // skip spaces and ':'
    while (*str != '\0' && (isspace((unsigned char)*str) || *str == ':'))
    {
        ++str;
    }

    // collect token
    while (*str != '\0' && *str != ':' && !isspace((unsigned char)*str))
    {
        token += *str;
        ++str;
    }

    // 跳过后续的分隔符
    while (*str != '\0' && (isspace((unsigned char)*str) || *str == ':'))
        ++str;

    return token;
}

SCIP_RETCODE ReaderTSP::getNodesFromFile(
    SCIP_FILE *file,
    double *x_coords,
    double *y_coords,
    Graph *graph)
{
    char str[SCIP_MAXSTRLEN];
    int i = 0;
    int nodenumber;
    GRAPHNODE *node = &(graph->nodes[0]);

    while (i < graph->nnodes && !SCIPfeof(file))
    {
        // 读取下一行内容
        (void)SCIPfgets(str, SCIP_MAXSTRLEN, file);
        char *s = str;

        if (!SCIPstrToIntValue(str, &nodenumber, &s))
        {
            SCIPerrorMessage("Could not read node number:\n%s\n", str);
            return SCIP_INVALIDDATA;
        }

        if (!SCIPstrToIntValue(str, &nodenumber, &s))
        {
            SCIPerrorMessage("Could not read node number:\n%s\n", str);
            return SCIP_INVALIDDATA;
        }

        if (!SCIPstrToRealValue(s, &x_coords[i], &s))
        {
            SCIPerrorMessage("Could not read x coordinate:\n%s\n", str);
            return SCIP_INVALIDDATA;
        }

        if (!SCIPstrToRealValue(s, &y_coords[i], &s))
        {
            SCIPerrorMessage("Could not read y coordinate:\n%s\n", str);
            return SCIP_INVALIDDATA;
        }

        node->id = i;
        if (nodenumber - 1 != i)
        {
            cout << "warning: nodenumber <" << nodenumber << "> does not match its index in node list <" << i + 1
                 << ">. Node will get number " << i + 1 << " when naming variables and constraints!" << endl;
            return SCIP_INVALIDDATA;
        }
        node->x = x_coords[i];
        node->y = y_coords[i];
        node->first_edge = NULL;
        node++;
        i++;
    }

    assert(i == graph->nnodes);

    return SCIP_OKAY;
}

bool ReaderTSP::checkValid(
    GRAPH *graph,
    const std::string &name,
    const std::string &type,
    const std::string &edgeweighttype,
    int nnodes)
{
    // if something seems to be strange, it will be reported, that file was not valid
    if (nnodes < 1)
    {
        cout << "parse error in file <" << name << "> dimension should be greater than 0" << endl;
        return false;
    }

    if (type != "TSP")
    {
        cout << "parse error in file <" << name << "> type should be TSP" << endl;
        return false;
    }

    if (!(edgeweighttype == "EUC_2D" || edgeweighttype == "MAX_2D" || edgeweighttype == "MAN_2D" || edgeweighttype == "GEO" || edgeweighttype == "ATT"))
    {
        cout << "parse error in file <" << name << "> unknown weight type, should be EUC_2D, MAX_2D, MAN_2D, ATT, or GEO" << endl;
        return false;
    }

    if (graph == NULL)
    {
        cout << "error while reading file <" << name << ">, graph is uninitialized. ";
        cout << "Probably NODE_COORD_SECTION is missing" << endl;
        return false;
    }

    return true;
}

SCIP_RETCODE ReaderTSP::addVarToEdges(
    SCIP* scip,
    GRAPHEDGE* edge,
    SCIP_VAR* var
    )
{
    assert(scip != NULL);
    assert(edge != NULL);
    assert(var != NULL);

    edge->var = var;
    SCIP_CALL(SCIPcaptureVar(scip, edge->var));

    edge->back->var = edge->var;
    SCIP_CALL(SCIPcaptureVar(scip, edge->back->var));
    return SCIP_OKAY;
}

SCIP_DECL_READERFREE(ReaderTSP::scip_free)
{
    return SCIP_OKAY;
}

SCIP_DECL_READERREAD(ReaderTSP::scip_read)
{
    SCIP_RETCODE retcode;

    GRAPH *graph = NULL;
    GRAPHNODE *node;
    GRAPHNODE *nodestart;
    GRAPHNODE *nodeend;
    GRAPHEDGE *edgeforw;
    GRAPHEDGE *edgebackw;
    GRAPHEDGE *edge;

    double *x_coords = NULL;
    double *y_coords = NULL;

#ifdef SCIP_DEBUG
    double **weights = NULL;
#endif

    double x;
    double y;

    int nnodes = 0;
    int nedges = 0;
    int i;
    int j;

    string name = "MY_OWN_LITTLE_TSP";
    string token;
    string type = "TSP";
    string edgeweighttype = "EUC_2D";

    retcode = SCIP_OKAY;
    *result = SCIP_DIDNOTRUN;

    // 打开文件
    SCIP_FILE *file = SCIPfopen(filename, "r");
    if (!file)
    {
        return SCIP_READERROR;
    }

    // 读取第一行
    char str[SCIP_MAXSTRLEN];
    (void)SCIPfgets(str, SCIP_MAXSTRLEN, file);

    // 获取第一行的token
    char *s = str;
    token = getToken(s);

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
            if (!SCIPstrToIntValue(s, &nnodes, &s))
            {
                SCIPerrorMessage("Could not read number of nodes:\n%s\n", s);
                return SCIP_INVALIDDATA;
            }
            nedges = nnodes * (nnodes - 1);
        }
        else if (token == "EDGE_WEIGHT_TYPE")
        {
            edgeweighttype = getToken(s);
        }
        else if (token == "NODE_COORD_SECTION" || token == "DISPLAY_DATA_SECTION")
        {
            if (nnodes < 1)
            {
                retcode = SCIP_READERROR;
                break;
            }
            else if (create_graph(nnodes, nedges, &graph))
            {
                assert(x_coords == NULL);
                assert(y_coords == NULL);

                x_coords = new double[nnodes];
                y_coords = new double[nnodes];

                SCIP_CALL(getNodesFromFile(file, x_coords, y_coords, graph));
            }
            else
            {
                retcode = SCIP_NOMEMORY;
                break;
            }
        }
        else if (token == "COMMENT:" || token == "COMMENT" || token == "DISPLAY_DATA_TYPE" || token == "DISPLAY_DATA_TYPE:")
        {
            // do nothing
        }
        else if (token == "EOF")
        {
            break;
        }
        else if (token == "")
        {
            // do notning
        }
        else
        {
            cout << "parse error in file <" << name << "> unknown keyword <" << token << ">" << endl;
            return SCIP_READERROR;
        }

        (void)SCIPfgets(str, SCIP_MAXSTRLEN, file);
        s = str;
        token = getToken(s);
    } // 读取文件结束

    SCIPfclose(file);

    // 检查数据是否有问题
    if (!checkValid(graph, name, type, edgeweighttype, nnodes))
    {
        retcode = SCIP_READERROR;
    }

    assert(graph != NULL);

    if (retcode == SCIP_OKAY) {
        auto edgeforw = &(graph->edges[0]);
        auto edgebackw = &(graph->edges[nedges / 2]);
#ifdef SCIP_DEBUG
        weights = new double *[nnodes];
        for (i = 0; i < nnodes; ++i)
            weights[i] = new double[nnodes];
#endif
        // 创建图的边
        for (int i = 0; i < nnodes; i++) {
            auto nodestart = &graph->nodes[i];
            for (int j = i + 1; j < nnodes; j++) {
                auto nodeend = &graph->nodes[j];

                edgeforw->adjac = nodeend;
                edgebackw->adjac = nodestart;
                edgeforw->back = edgebackw;
                edgebackw->back = edgeforw;

                auto x = x_coords[nodestart->id] - x_coords[nodeend->id];
                auto y = y_coords[nodestart->id] - y_coords[nodeend->id];
                // 计算节点之间的距离
                if( edgeweighttype == "EUC_2D") {
                    edgeforw->length = sqrt( x*x + y*y );
                } else if( edgeweighttype == "MAX_2D") {
                    edgeforw->length = MAX( ABS(x), ABS(y) );
                } else if( edgeweighttype == "MAN_2D") {
                edgeforw->length = ABS(x) + ABS(y);
                } else if( edgeweighttype == "ATT") {
                edgeforw->length = ceil( sqrt( (x*x+y*y)/10.0 ) );
                } else if( edgeweighttype == "GEO") {
                const double pi =  3.141592653589793;
                double rads[4];
                double coords[4];
                double degs[4];
                double mins[4];
                double euler[3];
                int k;

                coords[0] = x_coords[(*nodestart).id]; /*lint !e613*/
                coords[1] = y_coords[(*nodestart).id]; /*lint !e613*/
                coords[2] = x_coords[(*nodeend).id]; /*lint !e613*/
                coords[3] = y_coords[(*nodeend).id]; /*lint !e613*/

                for( k = 0; k < 4; k++ )
                {
                    degs[k] = coords[k] >= 0 ? floor(coords[k]) : ceil(coords[k]);
                    mins[k] = coords[k] - degs[k];
                    rads[k] = pi*(degs[k]+5.0*mins[k]/3.0)/180.0;
                }

                euler[0] = cos(rads[1]-rads[3]);
                euler[1] = cos(rads[0]-rads[2]);
                euler[2] = cos(rads[0]+rads[2]);
                edgeforw->length = floor(6378.388 * acos(0.5*((1.0+euler[0])*euler[1]-(1.0-euler[0])*euler[2]))+1.0);
                }

                // 将边长四舍五入为整数
                if (round_lengths_) {
                    edgeforw->length = NINT(edgeforw->length);
                }
                edgebackw->length = edgeforw->length;
    #ifdef SCIP_DEBUG
                weights[i][j] = edgeforw->length;
                weights[j][i] = edgebackw->length;
    #endif
                // 将其中一条半有向边插入到该节点的边链表中
                if (nodestart->first_edge == NULL) {
                    nodestart->first_edge = edgeforw;
                    nodestart->first_edge->next = NULL;
                } else {
                    edgeforw->next = nodestart->first_edge;
                    nodestart->first_edge =edgeforw;
                }

                // dito
                if (nodeend->first_edge == NULL) {
                nodeend->first_edge = edgebackw;
                nodeend->first_edge->next = NULL;
                } else {
                edgebackw->next = nodeend->first_edge;
                nodeend->first_edge = edgebackw;
                }

                edgeforw++;
                edgebackw++;
            }
        }
    }

    delete[] x_coords;
    delete[] y_coords;

    if (retcode != SCIP_OKAY) {
#ifdef SCIP_DEBUG
        if (weights != NULL) {
            for (int i = 0; i < nnodes; i++) {
                delete[] weights[i];
            }
            delete[] weights;
        }
#endif
        return retcode;
    }
#ifdef SCIP_DEBUG
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < nnodes; i++) {
        for (j = 0; j < nnodes; j++) {
            std::cout << weights[i][j] << std::endl;
            delete[] weights[i];
        }
    }
    delete[] weights;
#endif
    
    SCIP_CALL(SCIPcreateObjProb(scip, name.c_str(), new ProbDataTSP(graph), TRUE));

    // 创建变量，并将其添加到边中
    for (int i = 0; i < nedges / 2; i++) {
        SCIP_VAR* var;
        stringstream varname;
        auto edge = &graph->edges[i];

        // 变量以其所代表的边链接的两个节点命名
        varname << "x_e_" << edge->back->adjac->id+1 << "-" << edge->adjac->id+1;
        SCIP_CALL(SCIPcreateVar(scip, &var, varname.str().c_str(), 0.0, 1.0, edge->length,
                SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
        
        // 将变量添加到SCIP以及图中
        SCIP_CALL(SCIPaddVar(scip, var));
        SCIP_CALL(addVarToEdges(scip, edge, var));

        SCIP_CALL(SCIPreleaseVar(scip, &var));
    }

    // 添加所有节点的度约束
    if (nnodes >= 2) {
        auto node = &(graph->nodes[0]);
        for (int i = 0; i < nnodes; i++, node++) {
            SCIP_CONS* cons;
            stringstream consname;
            consname << "deg_con_v" << node->id + 1;

            SCIP_CALL(SCIPcreateConsLinear(scip, &cons, consname.str().c_str(), 0, NULL, NULL, 2.0, 2.0,
                    TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
            auto edge = node->first_edge;

            while (edge != NULL) {
                SCIP_CALL(SCIPaddCoefLinear(scip, cons, edge->var, 1.0));
                edge = edge->next;
            }

            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));
        }
    }

    release_graph(&graph);
    *result = SCIP_SUCCESS;

    return SCIP_OKAY;
}

SCIP_DECL_READERWRITE(ReaderTSP::scip_write)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}