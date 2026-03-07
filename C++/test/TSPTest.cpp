#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>

#include "../src/BranchAndBound/TSP.h"

TEST(TSPConstructorTest, ValidInput) {
    int N = 3;
    std::vector<std::vector<double>> adj = {
        {0.0, 1.0, 2.0},
        {1.0, 0.0, 3.0},
        {2.0, 3.0, 0.0}
    };

    TSP tsp(N, adj);

    EXPECT_EQ(tsp.GetNodeNum(), N); 
    EXPECT_EQ(tsp.GetDistanceMatrix(), adj);
}

TEST(TSPTest, SearchPathTest) {
    std::vector<std::vector<double>> adj = { {0, 10, 15, 20},
        {10, 0, 35, 25},
        {15, 35, 0, 30},
        {20, 25, 30, 0}
    };

    TSP tsp(4, adj);

    tsp.SearchPath();

    EXPECT_EQ(tsp.GetPath(), std::vector<int>({0, 1, 3, 2, 0}));
}