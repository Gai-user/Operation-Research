#include <gtest/gtest.h>

#include <limits>
#include <vector>

#include "../src/Simplex/SimplexPrimal/SimplexPrimal.h"

namespace {

using simplex::LP;
using simplex::SimplexResult;
using simplex::SimplexPrimalSolver;
using simplex::SimplexStatus;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kTol = 1e-6;

// 构造默认 x >= 0 的 LP
LP makeLP(const std::vector<double>& cost, const std::vector<std::vector<double>>& rows,
          const std::vector<double>& rowLower, const std::vector<double>& rowUpper) {
    LP lp;
    lp.cost = cost;
    lp.rows = rows;
    lp.rowLower = rowLower;
    lp.rowUpper = rowUpper;
    lp.varLower.assign(cost.size(), 0.0);
    lp.varUpper.assign(cost.size(), kInf);
    return lp;
}

// 教科书 LP: min -3x - 2y, s.t. x + y <= 4, x + 3y <= 6, x,y >= 0
// 顶点 (0,0)/(4,0)/(3,1)/(0,2) 对应 3x+2y = 0/12/11/4，最优在 (4, 0), 目标值 -12
TEST(SimplexTest, TextbookLP) {
    const LP lp = makeLP({-3.0, -2.0}, {{1.0, 1.0}, {1.0, 3.0}}, {-kInf, -kInf}, {4.0, 6.0});
    const SimplexResult result = SimplexPrimalSolver(lp).solve();

    EXPECT_EQ(result.status, SimplexStatus::kOptimal);
    ASSERT_TRUE(result.objective.has_value());
    EXPECT_NEAR(*result.objective, -12.0, kTol);
    ASSERT_TRUE(result.solution.has_value());
    EXPECT_NEAR((*result.solution)[0], 4.0, kTol);
    EXPECT_NEAR((*result.solution)[1], 0.0, kTol);
}

// 需要 Phase I: min 2x + y, s.t. x + y >= 3, x <= 2, x,y >= 0 -> 最优 3, (x,y) = (0, 3)
TEST(SimplexTest, PhaseOneGreaterEqual) {
    const LP lp = makeLP({2.0, 1.0}, {{1.0, 1.0}, {1.0, 0.0}}, {3.0, -kInf}, {kInf, 2.0});
    const SimplexResult result = SimplexPrimalSolver(lp).solve();

    EXPECT_EQ(result.status, SimplexStatus::kOptimal);
    ASSERT_TRUE(result.objective.has_value());
    EXPECT_NEAR(*result.objective, 3.0, kTol);
    ASSERT_TRUE(result.solution.has_value());
    EXPECT_NEAR((*result.solution)[0], 0.0, kTol);
    EXPECT_NEAR((*result.solution)[1], 3.0, kTol);
}

// 等式约束 + 人工变量: min x + y, s.t. x + y = 5, x - y = 1, x,y >= 0 -> x=3, y=2, obj=5
TEST(SimplexTest, EqualityConstraints) {
    const LP lp = makeLP({1.0, 1.0}, {{1.0, 1.0}, {1.0, -1.0}}, {5.0, 1.0}, {5.0, 1.0});
    const SimplexResult result = SimplexPrimalSolver(lp).solve();

    EXPECT_EQ(result.status, SimplexStatus::kOptimal);
    ASSERT_TRUE(result.objective.has_value());
    EXPECT_NEAR(*result.objective, 5.0, kTol);
    ASSERT_TRUE(result.solution.has_value());
    EXPECT_NEAR((*result.solution)[0], 3.0, kTol);
    EXPECT_NEAR((*result.solution)[1], 2.0, kTol);
}

// Beale 退化/循环例子(经典反循环测试):
// min -3/4 x1 + 150 x2 - 1/50 x3 + 6 x4
// s.t. 1/4 x1 - 60 x2 - 1/25 x3 + 9 x4 <= 0
//      1/2 x1 - 90 x2 - 1/50 x3 + 3 x4 <= 0
//      x3 <= 1
// 最优值 -1/20 = -0.05
TEST(SimplexTest, BealeDegenerateCycling) {
    const LP lp = makeLP(
        {-0.75, 150.0, -0.02, 6.0},
        {{0.25, -60.0, -0.04, 9.0}, {0.5, -90.0, -0.02, 3.0}, {0.0, 0.0, 1.0, 0.0}},
        {-kInf, -kInf, -kInf}, {0.0, 0.0, 1.0});
    const SimplexResult result = SimplexPrimalSolver(lp).solve();

    EXPECT_EQ(result.status, SimplexStatus::kOptimal);
    ASSERT_TRUE(result.objective.has_value());
    EXPECT_NEAR(*result.objective, -0.05, kTol);
}

// 无界: min -x, s.t. x - y <= 1, x,y >= 0
TEST(SimplexTest, Unbounded) {
    const LP lp = makeLP({-1.0, 0.0}, {{1.0, -1.0}}, {-kInf}, {1.0});
    const SimplexResult result = SimplexPrimalSolver(lp).solve();
    EXPECT_EQ(result.status, SimplexStatus::kUnbounded);
}

// 不可行: x <= 1 且 x >= 3
TEST(SimplexTest, Infeasible) {
    const LP lp = makeLP({1.0}, {{1.0}, {1.0}}, {-kInf, 3.0}, {1.0, kInf});
    const SimplexResult result = SimplexPrimalSolver(lp).solve();
    EXPECT_EQ(result.status, SimplexStatus::kInfeasible);
}

// 一般变量界: min -x, s.t. x + y >= 2, x <= 3, x,y >= 0 -> x 取上界 3, obj = -3
TEST(SimplexTest, VariableUpperBound) {
    LP lp = makeLP({-1.0, 0.0}, {{1.0, 1.0}}, {2.0}, {kInf});
    lp.varUpper[0] = 3.0;
    const SimplexResult result = SimplexPrimalSolver(lp).solve();

    EXPECT_EQ(result.status, SimplexStatus::kOptimal);
    ASSERT_TRUE(result.objective.has_value());
    EXPECT_NEAR(*result.objective, -3.0, kTol);
    ASSERT_TRUE(result.solution.has_value());
    EXPECT_NEAR((*result.solution)[0], 3.0, kTol);
}

// 自由变量拆分: min 2x - y, s.t. y = 3, x + y >= 0, x,y 均为自由变量
// y = 3, x >= -3, 最优 x = -3 -> obj = -6 - 3 = -9
TEST(SimplexTest, FreeVariableSplit) {
    LP lp = makeLP({2.0, -1.0}, {{0.0, 1.0}, {1.0, 1.0}}, {3.0, 0.0}, {3.0, kInf});
    lp.varLower[0] = -kInf;
    lp.varLower[1] = -kInf;
    const SimplexResult result = SimplexPrimalSolver(lp).solve();

    EXPECT_EQ(result.status, SimplexStatus::kOptimal);
    ASSERT_TRUE(result.objective.has_value());
    EXPECT_NEAR(*result.objective, -9.0, kTol);
    ASSERT_TRUE(result.solution.has_value());
    EXPECT_NEAR((*result.solution)[0], -3.0, kTol);
    EXPECT_NEAR((*result.solution)[1], 3.0, kTol);
}

}  // namespace
