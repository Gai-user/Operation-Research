#ifndef OPERATIONSEARCH_SIMPLEX_SIMPLEX_PRIMAL_H_
#define OPERATIONSEARCH_SIMPLEX_SIMPLEX_PRIMAL_H_

#include <optional>
#include <string>
#include <vector>

namespace simplex
{

    // 原始单纯形法 -> 求解最小化问题
    //    min c^T x
    //    s.t. lhs <= A x <= rhs
    //         lower <= x <= upper
    // 无穷界用 std::numeric_limits<double>::infinity() 表示

    struct LP
    {
        std::vector<double> cost;              // 目标函数
        std::vector<double> varLower;          // 变量下界
        std::vector<double> varUpper;          // 变量上界
        std::vector<double> rowLower;          // 约束左端项
        std::vector<double> rowUpper;          // 约束右端项
        std::vector<std::vector<double>> rows; // 约束系数
    };

    enum class SimplexStatus
    {
        kOptimal,
        kUnbounded,
        kInfeasible,
        kIterationLimit
    };

    struct SimplexResult
    {
        SimplexStatus status = SimplexStatus::kIterationLimit;
        std::optional<std::vector<double>> solution; // 原变量取值（不包含松弛/人工变量）
        std::optional<double> objective;
        int iterations = 0;
        std::vector<int> basicIndex; // 最终基
    };

    /**
     * 修正单纯形求解器 v1.0
     *
     * 使用两阶段法求解LP问题
     *
     * 每次迭代：计算 y = c_B^T * B^-1 ----> 确定入基变量 pricing
     * ----> Ftran(alpha = B^-1 * a_j) ----> 确定出基变量 ratio test ----> 更新B^-1
     *
     * 周期性 refactorize 控制数值误差累积。
     *
     * 在v1.0版本中，显示维护稠密 B^-1，入基使用 Dantzig 规则，迭代停滞切换 Bland 规则防循环。
     */
    class SimplexPrimalSolver
    {
    public:
        explicit SimplexPrimalSolver(const LP &lp, bool trace = false);
        ~SimplexPrimalSolver() {};

        SimplexResult solve();

    private:
        struct LoopOutcome
        {
            SimplexStatus status = SimplexStatus::kIterationLimit;
            double objective = 0.0;
        };
        /**
         * 将 LP 转换为标准形式：
         *  min c^T * x
         *  s.t. A * x = b
         *       x >= 0
         * 模型中包含松弛/剩余/人工变量
         */
        void standardize();

        /**
         * 在当前基上跑单纯形循环
         * skipArtificial = true 是人工变量不参与 pricing （Phase 2）
         */
        LoopOutcome runLoop(const std::vector<double> &cost, bool skipArtificial);

        /* 计算基变量的值 x_B = B^-1 * b */
        std::vector<double> computeXb() const;

        /* 计算约束的对偶变量的值 y^T = c_B^T * B^-1*/
        std::vector<double> computeDual(const std::vector<double> &cost) const;

        /* 计算所有列的检验数 d_j = c_j - y^T * a_j */
        std::vector<double> computeReducedCost(const std::vector<double> &cost, const std::vector<double> &y) const;

        /* 计算 非基变量列 B^-1 * col */
        std::vector<double> ftran(const std::vector<double> &col) const;

        /* 从检验数中选出进基变量，找不到返回 -1 */
        int price(const std::vector<double> &d, bool skipArtificial, bool bland) const;

        /* 选择出基变量；ratiosOut 填入各行候选比值（alpha <= 容差 处置 inf） */
        int ratioTest(const std::vector<double> &xb, const std::vector<double> &alpha, bool bland,
                      std::vector<double> &ratiosOut) const;

        /* 更新 B^-1 */
        void updateInverse(int pivotRow, const std::vector<double> &alpha);

        /* 对当前基用 Gauss-Jordan 重算 B^-1 */
        bool refactorize();

        /* 去除人工变量 */
            void driveOutArtificials();

        bool isBasic(int col) const;

        void fillResult(SimplexResult &result) const;

        /* enter / leaveRow 为 -1 时表示"最优判定" / "无出基行" */
        void traceIteration(int enter, int leaveRow, double objective, const std::vector<double> &xb,
                            const std::vector<double> &y, const std::vector<double> &d,
                            const std::vector<double> &alpha, const std::vector<double> &ratios,
                            bool bland) const;

        /* 求解结束后的汇总输出 */
        void traceResult(const SimplexResult &result) const;

        const LP &lp_;
        bool trace_;
        int numVars = 0; // 原变量个数
        // 原变量 j 的变换常数: kShiftAndBound/kShift -> x = x' + shift; kFlip -> x = shift - x'
        std::vector<double> varShift_;

        // 标准型数据 （列存储：a_[j] 为模型中第j列的系数，长度m_)
        int m_ = 0;
        int n_ = 0;
        std::vector<std::vector<double>> a_;
        std::vector<double> b_;
        std::vector<double> cost_; // Phase II 目标系数
        std::vector<bool> isArtificial_;
        std::vector<std::string> colName_;
        // 列对原变量的归属: j+1 表示该列贡献 +x_j, -(j+1) 表示贡献 -x_j, 0 表示松弛/人工列
        std::vector<int> colOwner_;
        std::vector<int> basis_;
        std::vector<std::vector<double>> binv_; // m_ x m_
        int iterations_ = 0;
    };

}

#endif // OPERATIONSEARCH_SIMPLEX_SIMPLEX_PRIMAL_H_