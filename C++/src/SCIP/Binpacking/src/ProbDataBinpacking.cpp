#include "ProbDataBinpacking.h"

#include <cassert>
#include <cstring>
#include <utility>

#include "scip/cons_setppc.h"
#include "scip/scip.h"

#include "PricerBinpacking.h"
#include "VarDataBinpacking.h"

#define EVENTHDLR_NAME "addedvar"
#define EVENTHDLR_DESC "event handler for catching added variables"

namespace binpacking
{

    /* ------------------------------------------------------------------ */
    /* "addedvar" 事件回调（objscip 不提供事件封装，保留静态 C 回调）        */
    /* ------------------------------------------------------------------ */

    /** SCIP_EVENTTYPE_VARADDED：把新变量加入问题数据 */
    static SCIP_DECL_EVENTEXEC(eventExecAddedVar)
    {
        assert(eventhdlr != nullptr);
        assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
        assert(event != nullptr);
        assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

        ProbDataBinpacking *probdata = static_cast<ProbDataBinpacking *>(SCIPgetObjProbData(scip));
        probdata->addVar(scip, SCIPeventGetVar(event));

        return SCIP_OKAY;
    }

    /* ------------------------------------------------------------------ */
    /* 构造 / 私有工具                                                      */
    /* ------------------------------------------------------------------ */

    ProbDataBinpacking::ProbDataBinpacking(std::vector<int> ids, std::vector<SCIP_Longint> weights, SCIP_Longint capacity,
                                           std::vector<SCIP_CONS *> conss)
        : conss_(std::move(conss)),
          ids_(std::move(ids)),
          weights_(std::move(weights)),
          nitems_(static_cast<int>(ids_.size())),
          capacity_(capacity)
    {
    }

    void ProbDataBinpacking::freeVarsConss(SCIP *scip)
    {
        for (SCIP_VAR *var : vars_)
        {
            SCIP_CALL_ABORT(SCIPreleaseVar(scip, &var));
        }

        for (SCIP_CONS *cons : conss_)
        {
            SCIP_CALL_ABORT(SCIPreleaseCons(scip, &cons));
        }

        vars_.clear();
        conss_.clear();
    }

    /** 初始列：每个物品单独装一个箱子的列（对照原 createInitialColumns）
     *
     * 原 C 顺序：createVar(先不设数据) → SCIPaddVar → SCIPprobdataAddVar → SCIPaddCoefSetppc
     * → 创建 vardata(consids=[i]) → SCIPvarSetData → SCIPchgVarUbLazy → SCIPreleaseVar。
     * C++ 中 vardata 在 SCIPcreateVarBinpacking 时一并传入（等价于"建列前就有数据"），其余顺序不变。
     */
    static SCIP_RETCODE createInitialColumns(SCIP *scip, ProbDataBinpacking *probdata)
    {
        const int nitems = probdata->nitems();
        const auto &ids = probdata->ids();
        const auto &conss = probdata->conss();

        for (int i = 0; i < nitems; ++i)
        {
            char name[SCIP_MAXSTRLEN];
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "item_%d", ids[i]);

            /* 该列只包含物品 i */
            SCIP_VAR *var = nullptr;
            SCIP_CALL(SCIPcreateVarBinpacking(scip, &var, name, 1.0, true, true, new VarDataBinpacking(std::vector<int>({i}))));

            SCIP_CALL(SCIPaddVar(scip, var));
            probdata->addVar(scip, var);
            SCIP_CALL(SCIPaddCoefSetppc(scip, conss[i], var));

            /* 上界设为 lazy：避免 LP 中 x<=1 产生可能为正的对偶变量 */
            SCIP_CALL(SCIPchgVarUbLazy(scip, var, 1.0));

            SCIP_CALL(SCIPreleaseVar(scip, &var));
        }

        return SCIP_OKAY;
    }

    /* ------------------------------------------------------------------ */
    /* objscip 回调                                                         */
    /* ------------------------------------------------------------------ */

    SCIP_RETCODE ProbDataBinpacking::scip_delorig(SCIP *scip)
    {
        freeVarsConss(scip);
        return SCIP_OKAY;
    }

    SCIP_RETCODE ProbDataBinpacking::scip_trans(SCIP *scip, scip::ObjProbData **objprobdata, SCIP_Bool *deleteobject)
    {
        /* 用相同数据构造变换后的问题数据（所有权交给 SCIP） */
        auto *transdata = new ProbDataBinpacking(ids_, weights_, capacity_);

        /* 拷贝原始约束/变量指针数组，再原位变换为 transformed 版本（对照原 probtransBinpacking） */
        transdata->conss_ = conss_;
        transdata->vars_ = vars_;
        SCIP_CALL(SCIPtransformConss(scip, transdata->nitems(), transdata->conss_.data(), transdata->conss_.data()));
        SCIP_CALL(SCIPtransformVars(scip, transdata->nvars(), transdata->vars_.data(), transdata->vars_.data()));

        *objprobdata = transdata;
        *deleteobject = TRUE;
        return SCIP_OKAY;
    }

    SCIP_RETCODE ProbDataBinpacking::scip_deltrans(SCIP *scip)
    {
        freeVarsConss(scip);
        return SCIP_OKAY;
    }

    SCIP_RETCODE ProbDataBinpacking::scip_initsol(SCIP *scip)
    {
        SCIP_EVENTHDLR *eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
        assert(eventhdlr != nullptr);
        SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, nullptr, nullptr));
        return SCIP_OKAY;
    }

    SCIP_RETCODE ProbDataBinpacking::scip_exitsol(SCIP *scip, SCIP_Bool /*restart*/)
    {
        SCIP_EVENTHDLR *eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
        assert(eventhdlr != nullptr);
        SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, nullptr, -1));
        return SCIP_OKAY;
    }

    /* ------------------------------------------------------------------ */
    /* 接口方法                                                             */
    /* ------------------------------------------------------------------ */
    SCIP_RETCODE SCIPprobdataCreate(
        SCIP *scip,
        const char *probname,
        const std::vector<int> &ids,
        const std::vector<SCIP_Longint> &weights,
        SCIP_Longint capacity)
    {
        const int nitems = static_cast<int>(ids.size());

        /* 若无 "addedvar" 事件处理器则创建 */
        if (SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == nullptr)
            SCIP_CALL(SCIPincludeEventhdlrBasic(scip, nullptr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, nullptr));

        /* 先创建问题再创建约束：SCIPcreateConsSetcover 在 PROBLEM 阶段才走"原始问题"分支，
         * 若在更早阶段调用会被误判为"变换后问题"而要求变换变量（对照原 SCIPcreateProbBasic 的先后顺序）。
         * 问题对象在创建时以空约束数组注册，随后通过 setConss 挂入。 */
        auto* probdata = new ProbDataBinpacking(ids, weights, capacity);
        SCIP_CALL(SCIPcreateObjProb(scip, probname, probdata, true));

        /* 设置目标函数 */
        SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));
        SCIP_CALL(SCIPsetObjIntegral(scip));

        /* 为每个物品创建 set covering 约束 （modifiable=true, 允许列生成加入新列）*/
        std::vector<SCIP_CONS*> conss(nitems, nullptr);
        for (int i = 0; i < nitems; ++i)
        {
            char name[SCIP_MAXSTRLEN];
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "item_%d", ids[i]);
            SCIP_CALL(SCIPcreateConsBasicSetcover(scip, &conss[i], name, 0, nullptr));
            SCIP_CALL(SCIPsetConsModifiable(scip, conss[i], true));
            SCIP_CALL(SCIPaddCons(scip, conss[i]));
        }

        probdata->setConss(std::move(conss));

        /* 创建初始列，将变量挂载到 probdata*/
        SCIP_CALL(createInitialColumns(scip, probdata));

        /* 激活定价器 - 创建子问题并求解 */
        SCIP_CALL(SCIPpricerBinpackingActivate(scip, probdata->conss(), weights, ids, capacity));

        return SCIP_OKAY;
    }

    void ProbDataBinpacking::addVar(SCIP* scip, SCIP_VAR* var)
    {
        SCIP_CALL_ABORT(SCIPcaptureVar(scip, var));
        vars_.push_back(var);
    }

} // namespace binpacking