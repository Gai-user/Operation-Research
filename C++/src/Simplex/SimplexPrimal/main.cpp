#include <limits>
#include <vector>

#include "SimplexPrimal.h"

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

int main() {
    const LP lp = makeLP({2.0, 1.0}, {{1.0, 1.0}, {1.0, 0.0}}, {3.0, -kInf}, {kInf, 2.0});
    const SimplexResult result = SimplexPrimalSolver(lp, true).solve();
    return 0;
}