#pragma once

#include <vector>

#include "objscip/objscip.h"

namespace binpacking
{

    /**
     * @brief 装箱问题列生成定价器
     *
     * 定价子问题 -> 0-1背包问题：在容量约束下最大化物品对偶值之和
     */
    class PricerBinpacking : public scip::ObjPricer
    {
    public:
        PricerBinpacking(
            SCIP *scip,
            const std::vector<SCIP_CONS *> &conss,
            const std::vector<SCIP_Longint> &weights,
            const std::vector<int> &ids,
            SCIP_Longint capacity);

        ~PricerBinpacking() override = default;

        virtual SCIP_DECL_PRICERFREE(scip_free) override;
        virtual SCIP_DECL_PRICERINIT(scip_init) override;
        virtual SCIP_DECL_PRICEREXITSOL(scip_exitsol) override;
        virtual SCIP_DECL_PRICERREDCOST(scip_redcost) override;
        virtual SCIP_DECL_PRICERFARKAS(scip_farkas) override;

    private:
        SCIP_CONSHDLR *conshdlr_;          // "samediff" 约束处理器
        std::vector<SCIP_CONS *> conss_;   // 每个物品一个set covering 约束
        std::vector<SCIP_Longint> weights_; // 物品重量
        std::vector<int> ids_;              // 物品 ID
        int nitems_;                        // 物品数量
        SCIP_Longint capacity_;            // 箱子容量

        SCIP_RETCODE doPricing(SCIP *scip, SCIP_Bool farkas, SCIP_RESULT *result);
        SCIP_RETCODE initPricing(SCIP *scip, SCIP_Bool isfrakas, SCIP *subscip, SCIP_VAR **vars);
        SCIP_RETCODE addBranchingDecisionConss(SCIP *scip, SCIP *subscip, SCIP_VAR **vars);
        SCIP_RETCODE addFixedVarsConss(SCIP *scip, SCIP *subscip, SCIP_VAR **vars);
    };

    /** 将问题特定数据交给定价器并激活（原 SCIPpricerBinpackingActivate）
     *
     * 构造 PricerBinpacking 对象，SCIPincludeObjPricer(deleteobject=TRUE) 后激活。
     */
    SCIP_RETCODE SCIPpricerBinpackingActivate(
        SCIP *scip,
        const std::vector<SCIP_CONS *> &conss,
        const std::vector<SCIP_Longint> &weights,
        const std::vector<int> &ids,
        SCIP_Longint capacity);

} // namespace binpacking