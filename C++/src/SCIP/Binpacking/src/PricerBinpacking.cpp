#include "PricerBinpacking.h"

#include <cassert>
#include <cstring>
#include <string>
#include <vector>

#include "scip/cons_knapsack.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "ConshdlrSamediff.h"
#include "ProbDataBinpacking.h"
#include "VarDataBinpacking.h"

#define PRICER_NAME "binpacking"
#define PRICER_DESC "pricer for binpacking tours"
#define PRICER_PRIORITY 0
#define PRICER_DELAY TRUE /* only call pricer if all problem variables have non-negative reduced costs */

namespace binpacking
{

    /* 构造 */
    PricerBinpacking::PricerBinpacking(
        SCIP *scip,
        const std::vector<SCIP_CONS *> &conss,
        const std::vector<SCIP_Longint> &weights,
        const std::vector<int> &ids,
        SCIP_Longint capacity) : ObjPricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY),
                                 conshdlr_(SCIPfindConshdlr(scip, "samediff")),
                                 conss_(conss),
                                 weights_(weights),
                                 ids_(ids),
                                 nitems_(static_cast<int>(conss.size())),
                                 capacity_(capacity)
    {
        assert(conshdlr_ != nullptr);
        assert(nitems_ > 0);
    }

    /** 把当前节点路径上的 SAME/DIFFER 分支决策转成定价子问题的约束
     *
     * 对每个活动的分支决策约束：
     *   - SAME：x1 = x2  ⇔  0 <= x1 - x2 <= 0，即 varbound x1 - x2 <= 0 且 x1 - x2 >= 0（lhs=rhs=0, coef=-1）
     *   - DIFFER：x1 + x2 <= 1  ⇔  -inf <= x1 + x2 <= 1（lhs=-inf, rhs=1, coef=+1）
     */
    SCIP_RETCODE PricerBinpacking::addBranchingDecisionConss(
        SCIP* scip,
        SCIP* subscip,
        SCIP_VAR** vars)
    {
        assert(scip != nullptr);
        assert(subscip != nullptr);
        assert(conshdlr_ != nullptr);

        /* 收集所有分支决策约束 */
        SCIP_CONS** conss = SCIPconshdlrGetConss(conshdlr_);
        const int nconss = SCIPconshdlrGetNConss(conshdlr_);

        for (int c = 0; c < nconss; ++c)
        {
            SCIP_CONS* cons = conss[c];
            
            /* 若是当前节点该约束为激活，跳过 */
            if (!SCIPconsIsActive(cons)) { continue; }

            const int id1 = SCIPgetItemid1Samediff(cons);
            const int id2 = SCIPgetItemid2Samediff(cons);
            const ConsType type = SCIPgetTypeSamediff(cons);

            /* 根据分支类型，添加分支约束 */
            SCIP_Real vbdcoef;
            SCIP_Real lhs;
            SCIP_Real rhs;

            if (type == ConsType::Same)
            {
                lhs = 0.0;
                rhs = 0.0;
                vbdcoef = -1.0;
            }
            else if (type == ConsType::Differ)
            {
                lhs = -SCIPinfinity(scip);
                rhs = 1.0;
                vbdcoef = 1.0;
            }
            else
            {
                SCIPerrorMessage("unknown constraint type <%d>\n", static_cast<int>(type));
                return SCIP_INVALIDDATA;
            }

            SCIP_CONS* varbound = nullptr;
            SCIP_CALL(SCIPcreateConsBasicVarbound(subscip, &varbound, SCIPconsGetName(cons),
                                                 vars[id1], vars[id2], vbdcoef, lhs, rhs));
            SCIP_CALL(SCIPaddCons(subscip, varbound));
            SCIP_CALL(SCIPreleaseCons(subscip, &varbound));
        }

        return SCIP_OKAY;
    }

    /** 避免重复生成被局部固定为 0 的列
     *
     * 对被固定为 0 的装箱变量，其对应的物品组合不该再次被定价子问题生成。对每个这样的变量，
     * 构造一个 logicor 约束，其中包含该组合的所有物品变量（选中该组合的物品即变量=1，故对
     * 物品集合 S 取逻辑或 ¬x_i 禁止 S 全选）。
     */
    SCIP_RETCODE PricerBinpacking::addFixedVarsConss(
        SCIP *scip,
        SCIP *subscip,
        SCIP_VAR **vars)
    {
        SCIP_VAR **origvars = SCIPgetVars(scip);
        const int norigvars = SCIPgetNVars(scip);

        for (int v = 0; v < norigvars; ++v)
        {
            assert(SCIPvarGetType(origvars[v]) == SCIP_VARTYPE_BINARY);

            /* 上界 < 0.5 的 0-1 变量，即认为固定位0 */
            if (SCIPvarGetUbLocal(origvars[v]) < 0.5)
            {
                VarDataBinpacking *vardata = static_cast<VarDataBinpacking *>(SCIPgetObjVardata(scip, origvars[v]));
                const std::vector<int> &consids = vardata->consids();
                const int nconsids = vardata->nconsids();
                SCIP_Bool needed = true;

                std::vector<SCIP_VAR *> logicorvars(nitems_, nullptr);
                int nlogicorvars = 0;
                int c = 0;
                int consid = consids[0];
                int nvars = 0;

                /* 遍历所有物品约束；若该组合中的某个物品约束已禁用（物品已固定装进某箱），则无需再加约束 */
                for (int o = 0; o < nitems_ && needed; ++o)
                {
                    assert(o <= consid);
                    SCIP_CONS *cons = conss_[o];
                    if (SCIPconsIsEnabled(cons))
                    {
                        assert(SCIPgetNFixedonesSetppc(scip, cons) == 0);
                        SCIP_VAR *var = vars[nvars];
                        ++nvars;
                        assert(var != nullptr);
                        /* 该物品属于被禁止的组合，取反*/
                        if (o == consid)
                        {
                            SCIP_CALL(SCIPgetNegatedVar(subscip, var, &var));
                        }

                        logicorvars[nlogicorvars] = var;
                        ++nlogicorvars;
                    }
                    else if (o == consid)
                    {
                        needed = FALSE;
                    }

                    if (o == consid)
                    {
                        ++c;
                        if (c == nconsids)
                        {
                            consid = nitems_ + 100;
                        }
                        else
                        {
                            assert(consid < consids[c]);
                            consid = consids[c];
                        }
                    }
                }

                if (needed)
                {
                    SCIP_CONS *cons = nullptr;
                    SCIP_CALL(SCIPcreateConsBasicLogicor(subscip, &cons, SCIPvarGetName(origvars[v]), nlogicorvars, logicorvars.data()));
                    SCIP_CALL(SCIPsetConsInitial(subscip, cons, FALSE));
                    SCIP_CALL(SCIPaddCons(subscip, cons));
                    SCIP_CALL(SCIPreleaseCons(subscip, &cons));
                }
            }
        }

        return SCIP_OKAY;
    }
    /**
     * 初始化定价子问题，创建变量，添加约束
     * 物品已装进某个箱子或被固定位1， 跳过
     */
    SCIP_RETCODE PricerBinpacking::initPricing(
        SCIP *scip,
        SCIP_Bool isfarkas,
        SCIP *subscip,
        SCIP_VAR **vars)
    {
        assert(SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM);

        std::vector<SCIP_Longint> vals(nitems_, 0);
        int nvars = 0;

        for (int c = 0; c < nitems_; ++c)
        {
            SCIP_CONS *cons = conss_[c];
            assert(strncmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "setppc", 6) == 0);

            /* 跳过禁用的约束 */
            if (!SCIPconsIsEnabled(cons))
            {
                continue;
            }

            if (SCIPgetNFixedonesSetppc(scip, cons) == 1)
            {
                SCIP_CALL(SCIPdelConsLocal(scip, cons));
                continue;
            }

            /* 获取主问题对偶值 */
            const SCIP_Real dual = isfarkas ? SCIPgetDualfarkasSetppc(scip, cons) : SCIPgetDualsolSetppc(scip, cons);

            SCIP_VAR *var = nullptr;
            SCIP_CALL(SCIPcreateVarBasic(subscip, &var, SCIPconsGetName(cons), 0.0, 1.0, dual, SCIP_VARTYPE_BINARY));
            SCIP_CALL(SCIPaddVar(subscip, var));

            vals[nvars] = weights_[c];
            vars[nvars] = var;
            ++nvars;

            SCIP_CALL(SCIPreleaseVar(subscip, &var));
        }

        /* 容量约束 */
        SCIP_CONS *cons = nullptr;
        SCIP_CALL(SCIPcreateConsBasicKnapsack(subscip, &cons, "capacity", nvars, vars, vals.data(), capacity_));
        SCIP_CALL(SCIPaddCons(subscip, cons));
        SCIP_CALL(SCIPreleaseCons(subscip, &cons));

        /* 分支决策约束 + 避免重生成固定为0的列*/
        SCIP_CALL(addBranchingDecisionConss(scip, subscip, vars));
        SCIP_CALL(addFixedVarsConss(scip, subscip, vars));

        return SCIP_OKAY;
    }

    /* 执行列生成：求解一个sub SCIP -> 0-1背包问题 */
    SCIP_RETCODE PricerBinpacking::doPricing(SCIP *scip, SCIP_Bool isfarkas, SCIP_RESULT *result)
    {
        assert(scip != nullptr);
        assert(result != nullptr);

        *result = SCIP_DIDNOTRUN;

        /* 剩余时间/内存限制传给 sun scip */
        SCIP_Real timelimit;
        SCIP_CALL(SCIPgetRealParam(scip, "limits/time", &timelimit));
        if (!SCIPisInfinity(scip, timelimit))
        {
            timelimit -= SCIPgetSolvingTime(scip);
        }
        SCIP_Real memorylimit;
        SCIP_CALL(SCIPgetRealParam(scip, "limits/memory", &memorylimit));
        if (!SCIPisInfinity(scip, memorylimit))
        {
            memorylimit -= SCIPgetMemUsed(scip) / 1048576.0;
        }

        /* 创建并配置 sub SCIP */
        SCIP *subscip = nullptr;
        SCIP_CALL(SCIPcreate(&subscip));
        SCIP_CALL(SCIPincludeDefaultPlugins(subscip));
        SCIP_CALL(SCIPcreateProbBasic(subscip, "pricing"));
        SCIP_CALL(SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE));

        /* 子问题不响应 CTRL-C、关闭输出、设置限制 */
        SCIP_CALL(SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE));
        SCIP_CALL(SCIPsetIntParam(subscip, "display/verblevel", 0));
        SCIP_CALL(SCIPsetRealParam(subscip, "limits/time", timelimit));
        SCIP_CALL(SCIPsetRealParam(subscip, "limits/memory", memorylimit));

        /* 初始化定价子问题 - vars 只填前nvars个，其余为nullptr */
        std::vector<SCIP_VAR *> vars(nitems_, nullptr);
        SCIP_CALL(initPricing(scip, isfarkas, subscip, vars.data()));

        /* 求解 sun scip*/
        SCIP_CALL(SCIPsolve(subscip));

        const int nsols = SCIPgetNSols(subscip);
        SCIP_SOL **sols = SCIPgetSols(subscip);
        SCIP_Bool addvar = false;

        /* 遍历子问题的所有解：目标值 > 1.0 （redcost）或者 > 0.0 (Farkas)，则生成新列 */
        for (int s = 0; s < nsols; ++s)
        {
            SCIP_Bool feasible;
            assert(s == 0 || SCIPisFeasGE(subscip, SCIPgetSolOrigObj(subscip, sols[s - 1]), SCIPgetSolOrigObj(subscip, sols[s])));

            SCIP_CALL(SCIPcheckSolOrig(subscip, sols[s], &feasible, false, false));

            if (!feasible)
            {
                SCIPwarningMessage(scip, "solution in pricing problem (capacity <%" SCIP_LONGINT_FORMAT ">) is infeasible\n", capacity_);
                continue;
            }

            if (SCIPisFeasGT(subscip, SCIPgetSolOrigObj(subscip, sols[s]), isfarkas ? 0.0 : 1.0))
            {
                /* 收集新列中包含的物品 */
                std::vector<int> consids;
                consids.reserve(nitems_);
                std::string name = "items";

                for (int o = 0, v = 0; o < nitems_; ++o)
                {
                    if (!SCIPconsIsEnabled(conss_[o]))
                    {
                        continue;
                    }

                    assert(SCIPgetNFixedonesSetppc(scip, conss_[o]) == 0);

                    if (SCIPgetSolVal(subscip, sols[s], vars[v]) > 0.5)
                    {
                        name += "_" + std::to_string(ids_[o]);
                        consids.push_back(o);
                    }
                    else
                    {
                        assert(SCIPisFeasEQ(subscip, SCIPgetSolVal(subscip, sols[s], vars[v]), 0.0));
                    }

                    ++v;
                }
                /* 创建新列 */
                SCIP_VAR *var = nullptr;
                SCIP_CALL(SCIPcreateVarBinpacking(scip, &var, name.c_str(), 1.0, false, true, new VarDataBinpacking(consids)));

                /* 将新列添加到主问题中 */
                SCIP_CALL(SCIPaddPricedVar(scip, var, 1.0));
                addvar = true;

                /* 上界设为 lazy，避免 LP 中 x<=1 产生可能为正的对偶变量 */
                SCIP_CALL(SCIPchgVarUbLazy(scip, var, 1.0));

                /* 更新约束中的列 */
                for (int o : consids)
                {
                    assert(SCIPconsIsEnabled(conss_[o]));
                    SCIP_CALL(SCIPaddCoefSetppc(scip, conss_[o], var));
                }

                SCIP_CALL(SCIPreleaseVar(scip, &var));
            }
            else
            {
                break;
            }
        }

        if (addvar || SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL)
        {
            *result = SCIP_SUCCESS;
        }

        /* 释放sub scip*/
        SCIP_CALL(SCIPfree(&subscip));

        return SCIP_OKAY;
    }

    /* 回调实现 */
    SCIP_DECL_PRICERFREE(PricerBinpacking::scip_free)
    {
        /* conss_/weights_/ids_ 由 std::vector 自动回收；被 capture 的约束在 scip_exitsol 中 release */
        return SCIP_OKAY;
    }

    SCIP_DECL_PRICERINIT(PricerBinpacking::scip_init)
    {
        /* 把原始约束原位变换为 transformed 版本：release 原约束，取得变换后约束并 capture */
        for (int c = 0; c < nitems_; ++c)
        {
            SCIP_CONS *cons = conss_[c];
            SCIP_CALL(SCIPreleaseCons(scip, &conss_[c]));
            SCIP_CALL(SCIPgetTransformedCons(scip, cons, &conss_[c]));
            SCIP_CALL(SCIPcaptureCons(scip, conss_[c]));
        }

        return SCIP_OKAY;
    }

    SCIP_DECL_PRICEREXITSOL(PricerBinpacking::scip_exitsol)
    {
        /* 释放变换后 capture 的约束 */
        for (int c = 0; c < nitems_; ++c)
        {
            SCIP_CALL(SCIPreleaseCons(scip, &conss_[c]));
        }

        return SCIP_OKAY;
    }

    SCIP_DECL_PRICERREDCOST(PricerBinpacking::scip_redcost)
    {
        SCIP_CALL(doPricing(scip, FALSE, result));
        return SCIP_OKAY;
    }

    SCIP_DECL_PRICERFARKAS(PricerBinpacking::scip_farkas)
    {
        SCIP_CALL(doPricing(scip, TRUE, result));
        return SCIP_OKAY;
    }

    /* 接口方法 */
    SCIP_RETCODE SCIPpricerBinpackingActivate(
        SCIP *scip,
        const std::vector<SCIP_CONS *> &conss,
        const std::vector<SCIP_Longint> &weights,
        const std::vector<int> &ids,
        SCIP_Longint capacity)
    {
        assert(scip != nullptr);
        assert(!conss.empty());
        assert(!weights.empty());

        /* 为定价器 capture 所有 set covering 约束（所有权归定价器，scip_exitsol 中 release） */
        for (SCIP_CONS *cons : conss)
        {
            SCIP_CALL(SCIPcaptureCons(scip, cons));
        }

        /* 构造并包含定价器（deleteobject=TRUE：SCIP 释放 pricer 时销毁对象） */
        SCIP_CALL(SCIPincludeObjPricer(scip, new PricerBinpacking(scip, conss, weights, ids, capacity), true));

        /* 激活定价器 */
        SCIP_CALL(SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME)));

        return SCIP_OKAY;
    }

} // namespace binpacking