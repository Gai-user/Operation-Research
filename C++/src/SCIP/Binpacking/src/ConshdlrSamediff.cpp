#include "ConshdlrSamediff.h"

#include <algorithm>
#include <cassert>
#include <cstring>

#include "scip/scip.h"

#include "ProbDataBinpacking.h"
#include "VarDataBinpacking.h"

#define CONSHDLR_NAME "samediff"
#define CONSHDLR_DESC "stores the local branching decisions"

namespace binpacking
{
    /** 若是某一装箱方案违反分支约束，对应的变量局部固定为0 */
    static SCIP_RETCODE checkVariable(
        SCIP *scip,
        SCIP_ConsData *consdata,
        SCIP_VAR *var,
        int *nfixedvars,
        SCIP_Bool *cutoff)
    {
        *cutoff = false;

        /* 若是变量已固定为0，跳过*/
        if (SCIPvarGetUbLocal(var) < 0.5)
        {
            return SCIP_OKAY;
        }

        VarDataBinpacking *vardata = static_cast<VarDataBinpacking *>(SCIPgetObjVardata(scip, var));
        const std::vector<int> &consids = vardata->consids();

        /* consids是升序排列的，使用二分查找判断变量包含的物品中是否包含分支约束中的物品 */
        const SCIP_Bool existid1 = std::binary_search(consids.begin(), consids.end(), consdata->itemid1);
        const SCIP_Bool existid2 = std::binary_search(consids.begin(), consids.end(), consdata->itemid2);

        if ((consdata->type == ConsType::Same && existid1 != existid2) || (consdata->type == ConsType::Differ && existid1 == existid2))
        {
            SCIP_Bool infeasible = false;
            SCIP_Bool fixed = false;
            SCIP_CALL(SCIPfixVar(scip, var, 0.0, &infeasible, &fixed));

            if (infeasible)
            {
                assert(SCIPvarGetLbLocal(var) > 0.5);
                *cutoff = true;
            }
            else
            {
                assert(fixed);
                ++(*nfixedvars);
            }
            return SCIP_OKAY;
        }

        return SCIP_OKAY;
    }

    /* 检查新生成的变量是否违反分支约束 */
    static SCIP_RETCODE consdataFixVariables(
        SCIP *scip,
        SCIP_ConsData *consdata,
        const std::vector<SCIP_VAR *> &vars,
        SCIP_RESULT *result)
    {
        int nfixedvars = 0;
        SCIP_Bool cutoff = false;

        for (int v = consdata->npropagations; v < static_cast<int>(vars.size()) && !cutoff; ++v)
        {
            SCIP_CALL(checkVariable(scip, consdata, vars[v], &nfixedvars, &cutoff));
        }

        if (cutoff)
        {
            *result = SCIP_CUTOFF;
        }
        else if (nfixedvars > 0)
        {
            *result = SCIP_REDUCEDDOM;
        }

        return SCIP_OKAY;
    }

    /* 构造函数 */
    ConshdlrSamediff::ConshdlrSamediff(SCIP *scip)
        : ObjConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                      0,        /* sepapriority */
                      0,        /* enfopriority */
                      -9999999, /* checkpriority */
                      0,        /* sepafreq */
                      1,        /* propfreq */
                      1,        /* eagerfreq */
                      0,        /* maxprerounds */
                      FALSE,    /* delaysepa */
                      FALSE,    /* delayprop */
                      TRUE,     /* needscons */
                      SCIP_PROPTIMING_BEFORELP,
                      SCIP_PRESOLTIMING_FAST)
    {
    }

    /** 回调实现 */
    SCIP_RETCODE ConshdlrSamediff::scip_delete(SCIP * /*scip*/, SCIP_CONSHDLR * /*conshdlr*/, SCIP_CONS * /*cons*/, SCIP_CONSDATA **consdata)
    {
        assert(consdata != nullptr);
        assert(*consdata != nullptr);

        delete *consdata;
        *consdata = nullptr;
        return SCIP_OKAY;
    }

    SCIP_RETCODE ConshdlrSamediff::scip_trans(
        SCIP *scip, SCIP_CONSHDLR *conshdlr, SCIP_CONS *sourcecons, SCIP_CONS **targetcons)
    {
        assert(scip != nullptr);
        assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
        assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
        assert(sourcecons != nullptr);
        assert(targetcons != nullptr);

        SCIP_ConsData *sourcedata = static_cast<SCIP_ConsData *>(SCIPconsGetData(sourcecons));
        assert(sourcedata != nullptr);

        /* 创建目标约束数据*/
        auto *targetdata = new SCIP_ConsData(sourcedata->itemid1, sourcedata->itemid2, sourcedata->type, sourcedata->node);

        /* 创建目标约束（复制全部 flag，对照 consTransSamediff） */
        SCIP_CALL(SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                                 SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                                 SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                                 SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
                                 SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)));

        return SCIP_OKAY;
    }

    SCIP_DECL_CONSENFOLP(ConshdlrSamediff::scip_enfolp)
    {
        /* 原 C 代码该回调为 NULL */
        return SCIP_OKAY;
    }

    SCIP_DECL_CONSENFOPS(ConshdlrSamediff::scip_enfops)
    {
        /* 原 C 代码该回调为 NULL */
        return SCIP_OKAY;
    }

    SCIP_DECL_CONSCHECK(ConshdlrSamediff::scip_check)
    {
        /* 原 C 代码该回调为 NULL */
        return SCIP_OKAY;
    }

    SCIP_RETCODE ConshdlrSamediff::scip_lock(
        SCIP * /*scip*/, SCIP_CONSHDLR * /*conshdlr*/, SCIP_CONS * /*cons*/, SCIP_LOCKTYPE /*locktype*/,
        int /*nlockspos*/, int /*nlocksneg*/)
    {
        /* 原 C 代码该回调为 NULL */
        return SCIP_OKAY;
    }

    SCIP_DECL_CONSPROP(ConshdlrSamediff::scip_prop)
    {
        assert(scip != nullptr);
        assert(result != nullptr);

        ProbDataBinpacking *probdata = static_cast<ProbDataBinpacking *>(SCIPgetObjProbData(scip));
        assert(probdata != nullptr);

        const auto &vars = probdata->vars();

        *result = SCIP_DIDNOTFIND;

        for (int c = 0; c < nconss; ++c)
        {
            SCIP_ConsData *consdata = static_cast<SCIP_ConsData *>(SCIPconsGetData(conss[c]));
            assert(consdata != nullptr);

            if (!consdata->propagated)
            {
                SCIP_CALL(consdataFixVariables(scip, consdata, vars, result));
                ++consdata->npropagations;

                if (*result != SCIP_CUTOFF)
                {
                    consdata->propagated = true;
                    consdata->npropagatedvars = static_cast<int>(vars.size());
                }
                else
                {
                    break;
                }
            }
        }
        return SCIP_OKAY;
    }

    SCIP_DECL_CONSACTIVE(ConshdlrSamediff::scip_active)
    {
        ProbDataBinpacking *probdata = static_cast<ProbDataBinpacking *>(SCIPgetObjProbData(scip));
        assert(probdata != nullptr);

        SCIP_ConsData *consdata = static_cast<SCIP_ConsData *>(SCIPconsGetData(cons));
        assert(consdata != nullptr);
        assert(consdata->npropagatedvars <= probdata->nvars());

        if (consdata->npropagatedvars != probdata->nvars())
        {
            consdata->propagated = false;
            SCIP_CALL(SCIPrepropagateNode(scip, consdata->node));
        }

        return SCIP_OKAY;
    }

    SCIP_DECL_CONSDEACTIVE(ConshdlrSamediff::scip_deactive)
    {
        /* 原 C 代码该回调为 NULL */
        return SCIP_OKAY;
    }

    SCIP_DECL_CONSPRINT(ConshdlrSamediff::scip_print)
    {
        /* 原 C 代码该回调为 NULL */
        return SCIP_OKAY;
    }

    SCIP_RETCODE SCIPcreateConsSamediff(
        SCIP *scip,
        SCIP_CONS **cons,
        const char *name,
        int itemid1,
        int itemid2,
        ConsType type,
        SCIP_NODE *node,
        SCIP_Bool local)
    {
        SCIP_CONSHDLR *conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

        if (conshdlr == nullptr)
        {
            SCIPerrorMessage("samediff constraint handler not found\n");
            return SCIP_PLUGINNOTFOUND;
        }

        auto *consdata = new SCIP_ConsData(itemid1, itemid2, type, node);

        SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
                                 local, FALSE, FALSE, FALSE, TRUE));

        return SCIP_OKAY;
    }

    int SCIPgetItemid1Samediff(SCIP_CONS *cons)
    {
        SCIP_ConsData *consdata = static_cast<SCIP_ConsData *>(SCIPconsGetData(cons));
        assert(consdata != nullptr);
        return consdata->itemid1;
    }

    int SCIPgetItemid2Samediff(SCIP_CONS *cons)
    {
        SCIP_ConsData *consdata = static_cast<SCIP_ConsData *>(SCIPconsGetData(cons));
        assert(consdata != nullptr);
        return consdata->itemid2;
    }

    ConsType SCIPgetTypeSamediff(SCIP_CONS *cons)
    {
        SCIP_ConsData *consdata = static_cast<SCIP_ConsData *>(SCIPconsGetData(cons));
        assert(consdata != nullptr);
        return consdata->type;
    }

} // namespace binpacking