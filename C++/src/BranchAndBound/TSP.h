/**
 * 使用分支定界求解TSP问题
 * 参考 https://www.geeksforgeeks.org/dsa/traveling-salesman-problem-using-branch-and-bound-2/
 */

#ifndef BRANCH_AND_BOUND_TSP_H
#define BRANCH_AND_BOUND_TSP_H
#include <vector>
#include <limits>

class TSP
{
public:
    TSP() = default;
    TSP(int N, const std::vector<std::vector<double>>& adj);
    ~TSP() = default;

    int GetNodeNum() const { return nodeNum; }
    std::vector<std::vector<double>> GetDistanceMatrix() { return distance; }
    void SearchPath();
    std::vector<int> GetPath() { return path;}

private:
    int nodeNum;
    const std::vector<std::vector<double>> distance; // 各点之间的距离
    std::vector<bool> visited;
    std::vector<int> path; // 旅行商的路径
    double finalResult = std::numeric_limits<double>::max(); // 最小权重

    bool isValidMatrix(int N, const std::vector<std::vector<double>>& adj);
    void CopyToFinal(const std::vector<int>& curPath);

    double FirstMinEdgeWeight(int i);
    double SecondMinEdgeWeight(int i);
    void BranchAndBound(double curBound, double curWeight, int level, std::vector<int>& curPath);

    
};

#endif //BRANCH_AND_BOUND_TSP_H