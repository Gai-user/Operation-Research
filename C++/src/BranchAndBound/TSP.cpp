#include <vector>
#include <stdexcept>
#include <limits>

#include "TSP.h"


TSP::TSP(int N, const std::vector<std::vector<double>>& adj) 
    : nodeNum(N), 
      visited(N, false), 
      distance(adj),
      path(N + 1)
{
    if (N < 0) {
        throw std::invalid_argument("The number of node must be positive!");
    }

    if (!isValidMatrix(N, adj)) {
        throw std::invalid_argument("Distance Matrix must be N x N");
    }
    
}

bool TSP::isValidMatrix(int N, const std::vector<std::vector<double>>& adj)
{
    if (adj.size() != N) {
        return false;
    }
    for (const auto& row : adj) {
        if (row.size() != N) {
            return false;
        }
    }
    return true;
}

void TSP::CopyToFinal(const std::vector<int>& curPath) 
{
    for (size_t i = 0; i < nodeNum; ++i) {
        path[i] = curPath[i];
    }
    path[nodeNum] = curPath[0];
}

/**
 * 寻找以i为端点的边的最小权重
 * @param i 端点编号
 * @return 最小权重
 */
double TSP::FirstMinEdgeWeight(int i)
{
    if (i < 0 || i >= nodeNum) {
        throw std::out_of_range("Invalid node index");
    }
    double minWeight = std::numeric_limits<double>::max();
    for (size_t j = 0; j < nodeNum; ++j) {
        if (distance[i][j] < minWeight && i != j) {
            minWeight = distance[i][j];
        }
    }
    return minWeight;
}

/**
 * 寻找以i为端点的边的第二小权重
 * @param i 端点编号
 * @return 次小权重
 */
double TSP::SecondMinEdgeWeight(int i)
{
    if (i < 0 || i >= nodeNum) {
        throw std::out_of_range("Invalid node index");
    }

    double firstMinWeight = std::numeric_limits<double>::max();
    double secondMinWeight = std::numeric_limits<double>::max();

    for (size_t j = 0; j < nodeNum; ++j) {
        if (i == j) {
            continue;
        }
        if (distance[i][j] < firstMinWeight) {
            secondMinWeight = firstMinWeight;
            firstMinWeight = distance[i][j];
        } else if (distance[i][j] < secondMinWeight && firstMinWeight != distance[i][j]) {
            secondMinWeight = distance[i][j];
        }
    }
    return secondMinWeight;
}

/**
 * 分支定界过程
 * @param curBound 当前节点的下界
 * @param curWeight 当前路径距离之和
 * @param level 在搜索树中所处的层数
 * @param curPath 当前解
 */
void TSP::BranchAndBound(double curBound, double curWeight, int level, std::vector<int>& curPath) 
{
    // 若是 level == N，表示所有节点都已访问
    if (level == nodeNum) {
        // 检查路径中最后一个顶点是否有一条边回到第一个顶点
        if (distance[curPath[level - 1]][curPath[0]] != 0) {
            double curResult = curWeight + distance[curPath[level - 1]][curPath[0]];
            if (curResult < finalResult) {
                finalResult = curResult;
                CopyToFinal(curPath);
            }
        }
        return;
    }

    // 对于其他层级，遍历所有顶点递归构建搜索树
    for (size_t i = 0; i < nodeNum; i++) {
        // 如果下一个顶点 i 不是当前顶点 (即排除 i, i)
        // 且尚未被访问过，则考虑选择它
        if (distance[curPath[level - 1]][i] != 0 && visited[i] == false) {
            double tempBound = curBound;
            curWeight += distance[curPath[level - 1]][i];
        
        
            // 这里计算下界的逻辑，没有想清楚
            if (level == 1) {
                curBound -= ((FirstMinEdgeWeight(curPath[level - 1]) + FirstMinEdgeWeight(i)) / 2);
            } else {
                curBound -= ((SecondMinEdgeWeight(curPath[level - 1]) + FirstMinEdgeWeight(i)) / 2);
            }

            // curBound + curWeight 是当前节点的实际下界
            // 下界小于已知的最优解，则继续探索该分支
            if (curBound + curWeight < finalResult) {
                curPath[level] = i;
                visited[i] = true;
                BranchAndBound(curBound, curWeight, level + 1, curPath);
                visited[i] = false;
            }
            
            // 否则，剪枝
            curWeight -= distance[curPath[level-1]][i];
            curBound = tempBound;
        }

    }

}

void TSP::SearchPath()
{
    double curBound = 0;
    std::vector<int> curPath(nodeNum, 0);
    visited.assign(nodeNum, false);

    // 计算初始下界
    for (size_t i = 0; i < nodeNum; i++) {
        curBound += (FirstMinEdgeWeight(i) + SecondMinEdgeWeight(i));
    }

    curBound /= 2;

    visited[0] = true;
    curPath[0] = 0;

    BranchAndBound(curBound, 0, 1, curPath);
}