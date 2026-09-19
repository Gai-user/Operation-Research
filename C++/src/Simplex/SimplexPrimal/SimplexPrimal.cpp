#include "SimplexPrimal.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string_view>
#include <utility>

namespace simplex
{

    namespace
    {
        constexpr double kPricingTol = 1e-7;  // 入基变量判负容差
        constexpr double kRatioTol = 1e-9;    // ratio test 主元容差
        constexpr double kZeroTol = 1e-7;     // 人工变量出基 / 可行性判定容差
        constexpr double kImproveTol = 1e-8;  // 目标值改善判定容差
        constexpr int kBlandThreshold = 20;   // 目标值停止多少次后切换Bland规则
        constexpr int kRefactorInterval = 50; // 重算 B^-1 的迭代次数
        constexpr int kMaxIterations = 10000;
        constexpr double kInf = std::numeric_limits<double>::infinity();

        double dot(const std::vector<double> &u, const std::vector<double> &v)
        {
            double s = 0.0;
            for (std::size_t i = 0; i < u.size(); ++i)
            {
                s += u[i] * v[i];
            }
            return s;
        }

        /* trace 用：定宽打印一个数值 */
        std::string fmt(const double v)
        {
            if (v == kInf || v == -kInf)
            {
                return (v > 0) ? "inf" : "-inf";
            }
            std::ostringstream os;
            os << std::fixed << std::setprecision(6) << v;
            return os.str();
        }

        /* trace 用：打印 [  1.000000   2.000000 ] 形式 */
        std::string formatVec(const std::vector<double> &v)
        {
            std::ostringstream os;
            os << "[";
            for (std::size_t i = 0; i < v.size(); ++i)
            {
                os << (i == 0 ? " " : "  ") << std::setw(10) << std::right << fmt(v[i]);
            }
            os << " ]";
            return os.str();
        }

        std::string_view statusName(const SimplexStatus status)
        {
            switch (status)
            {
            case SimplexStatus::kOptimal:
                return "optimal";
            case SimplexStatus::kUnbounded:
                return "unbounded";
            case SimplexStatus::kInfeasible:
                return "infeasible";
            case SimplexStatus::kIterationLimit:
                return "iteration limit";
            }
            return "unknown";
        }
    } // namespace

    SimplexPrimalSolver::SimplexPrimalSolver(const LP &lp, bool trace) : lp_(lp), trace_(trace) {}

    SimplexResult SimplexPrimalSolver::solve()
    {
        SimplexResult result;

        // 将模型转换为标准形式
        standardize();

        /* Phase 1：最小化人工变量之和*/
        if (std::any_of(isArtificial_.begin(), isArtificial_.end(), [](bool b)
                        { return b; }))
        {
            std::vector<double> phase1Cost(n_, 0.0);
            for (int j = 0; j < n_; ++j)
            {
                if (!isArtificial_[j])
                {
                    continue;
                }

                phase1Cost[j] = 1.0;
            }

            if (trace_)
            {
                std::cout << "======= Phase 1 =======" << std::endl;
            }

            // 单纯形迭代
            const LoopOutcome outcome = runLoop(phase1Cost, false);

            // 第一阶段求解没有找到最优解或者最优解大于0，原问题不可行
            if (outcome.status != SimplexStatus::kOptimal || outcome.objective > kZeroTol)
            {
                result.status = (outcome.status == SimplexStatus::kIterationLimit)
                                    ? SimplexStatus::kIterationLimit
                                    : SimplexStatus::kInfeasible;
                result.iterations = iterations_;
                if (trace_)
                {
                    traceResult(result);
                }
                return result;
            }

            driveOutArtificials();
        }

        /* ---------- Phase 2：根据第一阶段得到的基以及B^-1(没有人工变量)，求解问题*/
        if (trace_)
        {
            std::cout << "======= Phase 2 =======" << std::endl;
        }

        const LoopOutcome outcome = runLoop(cost_, true);

        result.status = outcome.status;
        result.iterations = iterations_;

        if (outcome.status == SimplexStatus::kOptimal)
        {
            fillResult(result);
        }

        if (trace_)
        {
            traceResult(result);
        }

        return result;
    }

    /**
     * 标准化模型：变量变换 + 行拆分 + 松弛/剩余/人工变量
     */
    void SimplexPrimalSolver::standardize()
    {
        numVars = static_cast<int>(lp_.cost.size());
        constexpr double kBndTol = 1e-12;

        // 变量变换表：记录原变量 j 的每种变换方式
        enum class VarKind
        {
            kShiftAndBound,
            kShift,
            kFlip,
            kFree
        };
        struct VarTransform
        {
            VarKind kind;
            double shift; // kShiftAndBound/kShift: x = x' + shift; kFlip: x = shift - x'
            double upper; // kShiftAndBound: x'的上界
        };
        std::vector<VarTransform> transform(numVars);

        // 新变量（原变量下标，符号）；自由变量拆分为两个；
        std::vector<std::pair<int, int>> varColumns;

        for (int j = 0; j < numVars; ++j)
        {
            const double lower = lp_.varLower[j];
            const double upper = lp_.varUpper[j];

            if (lower > -kInf && upper < kInf)
            {
                // l <= x <= u -> x = x' + l
                transform[j] = {VarKind::kShiftAndBound, lower, upper - lower};
                varColumns.emplace_back(j, 1);
            }
            else if (lower > -kInf)
            {
                // l <= x -> x = x' + l
                transform[j] = {VarKind::kShift, lower, kInf};
                varColumns.emplace_back(j, 1);
            }
            else if (upper < kInf)
            {
                // x <= u -> x = u - x'
                transform[j] = {VarKind::kFlip, upper, kInf};
                varColumns.emplace_back(j, -1);
            }
            else
            {
                // free -> x = x' - x''
                transform[j] = {VarKind::kFree, 0.0, kInf};
                varColumns.emplace_back(j, 1);
                varColumns.emplace_back(j, -1);
            }
        }

        // 记录各原变量的变换常数，供 fillResult 还原原变量取值
        varShift_.resize(numVars);
        for (int j = 0; j < numVars; ++j)
        {
            varShift_[j] = transform[j].shift;
        }

        // 约束变换 -> 将变量替换为更新后的变量，更新约束左端项与右端项
        struct ParsedRow
        {
            std::vector<double> coeff; // 新变量系数
            double lower;
            double upper;
        };
        const int numRows = static_cast<int>(lp_.rows.size());
        std::vector<ParsedRow> parsed;
        parsed.reserve(numRows);

        for (int i = 0; i < numRows; ++i)
        {
            ParsedRow row;
            row.coeff.assign(varColumns.size(), 0.0);
            row.lower = lp_.rowLower[i];
            row.upper = lp_.rowUpper[i];

            for (int slot = 0; slot < static_cast<int>(varColumns.size()); ++slot)
            {
                const auto [var, sign] = varColumns[slot];
                double c = lp_.rows[i][var];
                const VarTransform &tf = transform[var];

                switch (tf.kind)
                {
                case VarKind::kShiftAndBound:
                case VarKind::kShift:
                    // x = x' + shift: 常数项 c*shift 移动到 rhs，系数不变
                    row.lower -= c * tf.shift;
                    row.upper -= c * tf.shift;
                    break;
                case VarKind::kFlip:
                    // x = shift - x': 常数项 c*shift 移动到 rhs，系数取反
                    row.lower -= c * tf.shift;
                    row.upper -= c * tf.shift;
                    c = -c;
                    break;
                case VarKind::kFree:
                    c = (sign > 0) ? c : -c;
                    break;
                default:
                    break;
                }
                row.coeff[slot] = c;
            }
            parsed.push_back(std::move(row));
        }

        // 有界变量的上界：x' <= U，生成附加行
        for (int slot = 0; slot < static_cast<int>(varColumns.size()); ++slot)
        {
            const auto [var, sign] = varColumns[slot];
            const VarTransform &tf = transform[var];

            if (tf.kind != VarKind::kShiftAndBound)
            {
                continue;
            }

            ParsedRow row;
            row.coeff.assign(varColumns.size(), 0.0);
            row.coeff[slot] = 1.0;
            row.lower = -kInf;
            row.upper = tf.upper;
            parsed.push_back(std::move(row));
        }

        // 将约束拆分为单边约束，rhs >= 0
        enum class Sense
        {
            kLe,
            kGe,
            kEq
        };
        struct StdRow
        {
            std::vector<double> coeff;
            double rhs;
            Sense sense;
        };

        std::vector<StdRow> stdRows;

        for (ParsedRow &row : parsed)
        {
            const bool hasLo = row.lower > -kInf;
            const bool hasUp = row.upper < kInf;
            if (hasLo && hasUp)
            {
                if (row.upper - row.lower <= kBndTol)
                {
                    stdRows.push_back({row.coeff, row.lower, Sense::kEq});
                }
                else
                {
                    stdRows.push_back({row.coeff, row.lower, Sense::kGe});
                    stdRows.push_back({row.coeff, row.upper, Sense::kLe});
                }
            }
            else if (hasUp)
            {
                stdRows.push_back({row.coeff, row.upper, Sense::kLe});
            }
            else if (hasLo)
            {
                stdRows.push_back({row.coeff, row.lower, Sense::kGe});
            }

            // 两端均无穷的行，直接丢弃
        }

        for (StdRow &row : stdRows)
        {
            if (row.rhs < 0.0)
            {
                for (double &c : row.coeff)
                {
                    c = -c;
                }
                row.rhs = -row.rhs;

                if (row.sense == Sense::kLe)
                {
                    row.sense = Sense::kGe;
                }
                else if (row.sense == Sense::kGe)
                {
                    row.sense = Sense::kLe;
                }
            }
        }

        // 添加 松弛/剩余/人工变量
        m_ = static_cast<int>(stdRows.size());             // 约束数量
        n_ = static_cast<int>(varColumns.size()) + 2 * m_; // 变量数量 -> // 上界: GE 行最多加 2 列(剩余+人工), 预留后缩减

        a_.assign(n_, std::vector<double>(m_, 0.0)); // 变量在约束当中的系数，a[j] 表示第j个变量那一列
        b_.assign(m_, 0.0);                          // 约束右端项系数
        basis_.assign(m_, -1);                       // 基变量下标
        cost_.assign(n_, 0.0);                       // Phase II 目标系数
        colOwner_.assign(n_, 0);
        isArtificial_.assign(n_, false);                // 变量是否为人工变量
        binv_.assign(m_, std::vector<double>(m_, 0.0)); // B^-1 矩阵

        for (int i = 0; i < m_; ++i)
        {
            binv_[i][i] = 1.0;
        }

        colName_.assign(n_, "");

        for (int slot = 0; slot < static_cast<int>(varColumns.size()); ++slot)
        {
            const auto [var, sign] = varColumns[slot];
            for (int i = 0; i < m_; ++i)
            {
                a_[slot][i] = stdRows[i].coeff[slot];
            }
            cost_[slot] = lp_.cost[var] * sign;
            colOwner_[slot] = (var + 1) * sign;
            // 自由变量拆成两列，用 +/- 区分
            colName_[slot] = "x" + std::to_string(var + 1) +
                             ((transform[var].kind == VarKind::kFree) ? (sign > 0 ? "+" : "-") : "");
        }

        int col = static_cast<int>(varColumns.size());

        // 添加 松弛/剩余/人工变量
        for (int i = 0; i < m_; ++i)
        {
            const StdRow &row = stdRows[i];
            b_[i] = row.rhs;

            if (row.sense == Sense::kLe)
            {
                // 添加松弛变量：x <= rhs -> x + x' == rhs
                a_[col][i] = 1.0;
                basis_[i] = col;
                colName_[col] = "s" + std::to_string(i);
                ++col;
            }
            else if (row.sense == Sense::kGe)
            {
                // 添加剩余变量：x >= rhs -> x - x' == rhs
                a_[col][i] = -1.0;
                colName_[col] = "s" + std::to_string(i);
                ++col;
                // 添加人工变量，为了构造初始基矩阵：x - x' == rhs -> x - x' + x'' == rhs
                a_[col][i] = 1.0;
                isArtificial_[col] = true;
                basis_[i] = col;
                colName_[col] = "a" + std::to_string(i);
                ++col;
            }
            else
            {
                // 添加人工变量：x == rhs -> x + x' == rhs
                a_[col][i] = 1.0;
                isArtificial_[col] = true;
                basis_[i] = col;
                colName_[col] = "a" + std::to_string(i);
                ++col;
            }
        }

        n_ = col; // 实际列数
        a_.resize(n_);
        cost_.resize(n_);
        colOwner_.resize(n_);
        isArtificial_.resize(n_);
        colName_.resize(n_);
    }

    /* 主循环 */
    SimplexPrimalSolver::LoopOutcome SimplexPrimalSolver::runLoop(const std::vector<double> &cost, const bool skipArtificial)
    {
        LoopOutcome outcome;
        bool bland = false;
        double bestObj = kInf;
        int stagnant = 0;
        int phaseIters = 0;

        while (true)
        {
            if (++phaseIters > kMaxIterations)
            {
                outcome.status = SimplexStatus::kIterationLimit;
                return outcome;
            }
            ++iterations_;
            if (iterations_ % kRefactorInterval == 0)
            {
                refactorize();
            }

            // 计算基变量的值：B*x_b = b
            const std::vector<double> xb = computeXb();
            double obj = 0.0; // c_b * x_b
            for (int i = 0; i < m_; ++i)
            {
                obj += cost[basis_[i]] * xb[i];
            }

            outcome.objective = obj;

            if (obj < bestObj - kImproveTol)
            {
                bestObj = obj;
                stagnant = 0;
            }
            else if (++stagnant > kBlandThreshold)
            {
                bland = true; // 目标停滞，切 Bland 规则
            }

            // y = c_b * B^-1
            const std::vector<double> y = computeDual(cost);
            // 各非基变量的检验数
            const std::vector<double> reduced = computeReducedCost(cost, y);
            // 进基变量
            const int enter = price(reduced, skipArtificial, bland);

            if (enter < 0)
            {
                // 无负检验数 -> 当前基最优
                if (trace_)
                {
                    traceIteration(-1, -1, obj, xb, y, reduced, {}, {}, bland);
                }
                outcome.status = SimplexStatus::kOptimal;
                return outcome;
            }
            // 计算进基变量对应的列：new_a[j] = B^-1 * a[j]
            const std::vector<double> alpha = ftran(a_[enter]);
            // 获取最小比值 - 出基变量
            std::vector<double> ratios;
            const int leaveRow = ratioTest(xb, alpha, bland, ratios);
            if (leaveRow < 0)
            {
                if (trace_)
                {
                    traceIteration(enter, -1, obj, xb, y, reduced, alpha, ratios, bland);
                }
                outcome.status = SimplexStatus::kUnbounded;
                return outcome;
            }

            if (trace_)
            {
                traceIteration(enter, leaveRow, obj, xb, y, reduced, alpha, ratios, bland);
            }

            // 更新 B^-1
            updateInverse(leaveRow, alpha);
            basis_[leaveRow] = enter;
        }
    }

    /* 计算 基变量的 值*/
    std::vector<double> SimplexPrimalSolver::computeXb() const
    {
        // B * x_b = b - N * x_n -> x_b = b^-1 * b
        std::vector<double> xb(m_, 0.0);

        for (int i = 0; i < m_; ++i)
        {
            xb[i] = dot(binv_[i], b_);
        }
        return xb;
    }

    /* 计算 对偶乘子 */
    std::vector<double> SimplexPrimalSolver::computeDual(const std::vector<double> &cost) const
    {
        // y^T = c_B^T B^-1 => y_k = sum_i cost[basis_[i]] * binv_[i][k]
        std::vector<double> y(m_, 0.0);
        for (int i = 0; i < m_; ++i)
        {
            const double cb = cost[basis_[i]];
            if (cb == 0.0)
            {
                continue;
            }

            for (int k = 0; k < m_; ++k)
            {
                y[k] += cb * binv_[i][k];
            }
        }

        return y;
    }

    std::vector<double> SimplexPrimalSolver::computeReducedCost(const std::vector<double> &cost,
                                                               const std::vector<double> &y) const
    {
        // d_j = c_j - y^T * a_j
        std::vector<double> d(n_, 0.0);
        for (int j = 0; j < n_; ++j)
        {
            d[j] = cost[j] - dot(y, a_[j]);
        }
        return d;
    }

    std::vector<double> SimplexPrimalSolver::ftran(const std::vector<double> &col) const
    {
        std::vector<double> alpha(m_, 0.0);
        for (int i = 0; i < m_; ++i)
        {
            alpha[i] = dot(binv_[i], col);
        }
        return alpha;
    }

    bool SimplexPrimalSolver::isBasic(const int col) const
    {
        return std::find(basis_.begin(), basis_.end(), col) != basis_.end();
    }

    int SimplexPrimalSolver::price(const std::vector<double> &d, const bool skipArtificial,
                                   const bool bland) const
    {
        int best = -1;
        double bestReduced = -kPricingTol;
        for (int j = 0; j < n_; ++j)
        {
            if (skipArtificial && isArtificial_[j])
            {
                continue;
            }

            if (isBasic(j))
            {
                continue;
            }

            if (bland)
            {
                if (d[j] < -kPricingTol)
                {
                    return j;
                } // 最小下标
            }
            else if (d[j] < bestReduced)
            {
                bestReduced = d[j]; // 最负检验数
                best = j;
            }
        }
        return best;
    }

    int SimplexPrimalSolver::ratioTest(const std::vector<double> &xb, const std::vector<double> &alpha,
                                       const bool bland, std::vector<double> &ratiosOut) const
    {
        int bestRow = -1;
        double bestRatio = kInf;
        ratiosOut.assign(m_, kInf);
        for (int i = 0; i < m_; ++i)
        {
            if (alpha[i] <= kRatioTol)
            {
                continue;
            }
            const double ratio = xb[i] / alpha[i];
            ratiosOut[i] = ratio;

            if (ratio < bestRatio - kRatioTol)
            {
                bestRatio = ratio;
                bestRow = i;
            }
            else if (bland && ratio <= bestRatio + kRatioTol)
            {
                if (bestRow < 0 || basis_[i] < basis_[bestRow])
                {
                    bestRow = i;
                }
            }
        }
        return bestRow;
    }

    void SimplexPrimalSolver::updateInverse(const int pivotRow, const std::vector<double> &alpha)
    {
        const double piv = alpha[pivotRow];
        for (double &v : binv_[pivotRow])
        {
            v /= piv;
        }

        for (int i = 0; i < m_; ++i)
        {
            if (i == pivotRow)
            {
                continue;
            }

            const double factor = alpha[i];
            if (factor == 0.0)
            {
                continue;
            }
            for (int k = 0; k < m_; ++k)
            {
                binv_[i][k] -= factor * binv_[pivotRow][k];
            }
        }
    }

    bool SimplexPrimalSolver::refactorize()
    {
        // 增广 [B | I] 做 Gauss - Jordan，重新计算 B^-1，消除误差
        std::vector<std::vector<double>> aug(m_, std::vector<double>(2 * m_, 0.0));

        for (int col = 0; col < m_; ++col)
        {
            const std::vector<double> &acol = a_[basis_[col]];
            for (int i = 0; i < m_; ++i)
            {
                aug[i][col] = acol[i];
            }

            aug[col][m_ + col] = 1.0;
        }

        for (int col = 0; col < m_; ++col)
        {
            int piv = col;
            for (int r = col + 1; r < m_; ++r)
            {
                if (std::abs(aug[r][col]) > std::abs(aug[piv][col]))
                {
                    piv = r;
                }
            }

            if (std::abs(aug[piv][col]) < 1e-12)
            {
                return false;
            }

            std::swap(aug[col], aug[piv]);
            const double d = aug[col][col];
            for (int k = col; k < 2 * m_; ++k)
            {
                aug[col][k] /= d;
            }

            for (int r = 0; r < m_; ++r)
            {
                if (r == col)
                {
                    continue;
                }
                const double f = aug[r][col];
                for (int k = col; k < 2 * m_; ++k)
                {
                    aug[r][k] -= f * aug[col][k];
                }
            }
        }

        for (int i = 0; i < m_; ++i)
        {
            binv_[i].assign(aug[i].begin() + m_, aug[i].end());
        }

        return true;
    }

    /* Pahse 1 收尾：将仍留在基里的人工变量挤出； 挤不出的行是冗余行，删除 */
    void SimplexPrimalSolver::driveOutArtificials()
    {
        std::vector<int> redundantRows;
        for (int r = 0; r < m_; ++r)
        {
            if (!isArtificial_[basis_[r]])
            {
                continue;
            }
            bool pivoted = false;
            for (int j = 0; j < n_; ++j)
            {
                if (isArtificial_[j] || isBasic(j))
                {
                    continue;
                }
                const std::vector<double> &alpha = ftran(a_[j]);
                if (std::abs(alpha[r]) > kZeroTol)
                {
                    if (trace_)
                    {
                        std::cout << "[cleanup] 人工变量 " << colName_[basis_[r]] << " 出基, "
                                  << colName_[j] << " 入基\n";
                    }
                    updateInverse(r, alpha);
                    basis_[r] = j;
                    pivoted = true;
                    break; // 该行人工变量已出基，换下一行
                }
            }

            if (!pivoted)
            {
                redundantRows.push_back(r);
            }
        }

        if (redundantRows.empty())
        {
            return;
        }

        // 删除冗余行，重构标准型与B^-1
        std::vector<bool> drop(m_, false);
        for (const int r : redundantRows)
        {
            drop[r] = true;
        }
        std::vector<int> keep;
        for (int i = 0; i < m_; ++i)
        {
            if (!drop[i])
            {
                keep.push_back(i);
            }
        }

        for (auto &colData : a_)
        {
            std::vector<double> reduced;
            reduced.reserve(keep.size());
            for (const int i : keep)
            {
                reduced.push_back(colData[i]);
            }

            colData = std::move(reduced);
        }

        std::vector<double> newB;
        std::vector<int> newBasis;
        newB.reserve(keep.size());
        newBasis.reserve(keep.size());
        for (const int i : keep)
        {
            newB.push_back(b_[i]);
            newBasis.push_back(basis_[i]);
        }

        b_ = std::move(newB);
        basis_ = std::move(newBasis);
        m_ = static_cast<int>(keep.size());
        refactorize();
        if (trace_)
        {
            std::cout << "[cleanup] 删除 " << redundantRows.size() << " 个冗余约束行\n";
        }
    }

    void SimplexPrimalSolver::fillResult(SimplexResult &result) const
    {
        const std::vector<double> xb = computeXb();
        std::vector<double> xStd(n_, 0.0);
        for (int i = 0; i < m_; ++i)
            xStd[basis_[i]] = xb[i];

        std::vector<double> x(numVars, 0.0);
        for (int j = 0; j < n_; ++j)
        {
            if (colOwner_[j] == 0 || xStd[j] == 0.0)
                continue;
            const int owner = std::abs(colOwner_[j]) - 1;
            x[owner] += (colOwner_[j] > 0) ? xStd[j] : -xStd[j];
        }
        // 还原变换常数: x = x' + shift / x = shift - x'
        for (int j = 0; j < numVars; ++j)
        {
            x[j] += varShift_[j];
        }
        result.solution = std::move(x);
        result.objective = dot(lp_.cost, *result.solution);
        result.basicIndex = basis_;
    }

    void SimplexPrimalSolver::traceIteration(const int enter, const int leaveRow, const double objective,
                                            const std::vector<double> &xb, const std::vector<double> &y,
                                            const std::vector<double> &d, const std::vector<double> &alpha,
                                            const std::vector<double> &ratios, const bool bland) const
    {
        std::ostringstream xbLine;
        xbLine << "[";
        for (int i = 0; i < m_; ++i)
        {
            xbLine << (i == 0 ? " " : "  ") << colName_[basis_[i]] << "=" << fmt(xb[i]);
        }
        xbLine << " ]";

        std::cout << "[iter " << std::setw(4) << iterations_ << "] xB    = " << xbLine.str()
                  << "   obj = " << fmt(objective) << (bland ? "  [Bland]" : "") << '\n'
                  << "           y     = " << formatVec(y) << '\n'
                  << "           d     = " << formatVec(d) << '\n';

        if (enter < 0)
        {
            std::cout << "           enter = 无 (所有检验数 >= 0) -> optimal\n";
            return;
        }

        std::cout << "           enter = " << colName_[enter] << "  (d = " << fmt(d[enter]) << ")\n"
                  << "           alpha = " << formatVec(alpha) << '\n'
                  << "           ratio = " << formatVec(ratios) << '\n';

        if (leaveRow < 0)
        {
            std::cout << "           leave = 无 (alpha 无非正) -> " << colName_[enter]
                      << " 可无限增大 -> unbounded\n";
            return;
        }

        std::cout << "           leave = " << colName_[basis_[leaveRow]] << "  (ratio = "
                  << fmt(ratios[leaveRow]) << ")\n";
    }

    void SimplexPrimalSolver::traceResult(const SimplexResult &result) const
    {
        std::cout << "======= Result =======" << '\n'
                  << "status     = " << statusName(result.status) << '\n'
                  << "iterations = " << result.iterations << '\n';

        std::cout << "objective  = ";
        if (result.objective.has_value())
        {
            std::cout << fmt(*result.objective);
        }
        else
        {
            std::cout << "-";
        }
        std::cout << '\n';

        std::cout << "solution   = ";
        if (!result.solution.has_value())
        {
            std::cout << "-";
        }
        else
        {
            for (int j = 0; j < numVars; ++j)
            {
                if (j > 0)
                {
                    std::cout << ", ";
                }
                std::cout << "x" << (j + 1) << " = " << fmt((*result.solution)[j]);
            }
        }
        std::cout << '\n';
    }

} // simplex