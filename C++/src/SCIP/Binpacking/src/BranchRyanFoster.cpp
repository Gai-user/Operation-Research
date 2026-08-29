#include "BranchRyanFoster.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>

#include "scip/scip.h"

#include "ConshdlrSamediff.h"
#include "ProbDataBinpacking.h"
#include "VarDataBinpacking.h"

#define BRANCHRULE_NAME "RyanFoster"
#define BRANCHRULE_DESC "Ryan/Foster branching rule"
#define BRANCHRULE_PRIORITY 50000
#define BRANCHRULE_MAXDEPTH -1
#define BRANCHRULE_MAXBOUNDDIST 1.0

namespace binpacking
{

    BranchRyanFoster::BranchRyanFoster(SCIP *scip)
        : ObjBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                        BRANCHRULE_MAXBOUNDDIST)
    {
    }

    SCIP_DECL_BRANCHEXECLP(BranchRyanFoster::scip_execlp)
    {
        assert(scip != nullptr);
        assert(branchrule != nullptr);
        assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
        assert(result != nullptr);

        *result = SCIP_DIDNOTRUN;

        ProbDataBinpacking* probdata = static_cast<ProbDataBinpacking*>(SCIPgetObjProbData(scip));
        assert(probdata != nullptr);

        const int nitems = probdata->nitems();
        const auto& ids = probdata->ids();

        /* 三角矩阵：pairweights[i] 长度为 i+1，pairweights[i][j]（j<i）记"物品 j 与 i 同箱"的 LP 权重 */
        std::vector<std::vector<SCIP_Real>> pairweights(nitems);
        for (int i = 0; i < nitems; ++i)
        {
            pairweights[i].resize(i + 1, 0.0);
        }

        /* 获取当前 LP 的分支分支候选 */
        SCIP_VAR** lpcands = nullptr;
        SCIP_Real* lpcandsfrac = nullptr;
        int nlpcands = 0;
        SCIP_CALL(SCIPgetLPBranchCands(scip, &lpcands, nullptr, &lpcandsfrac, nullptr, &nlpcands, nullptr));
        assert(nlpcands > 0);

        /* 累加每个候选变量的 LP 值 */
        for (int v = 0; v < nlpcands; ++v)
        {
            const SCIP_Real solval = lpcandsfrac[v];

            VarDataBinpacking* vardata = static_cast<VarDataBinpacking*>(SCIPgetObjVardata(scip, lpcands[v]));
            const std::vector<int> consids = vardata->consids();
            const int nconsids = vardata->nconsids();
            assert(nconsids > 0);

            for (int i = 0; i < nconsids; ++i)
            {
                const int id1 = consids[i];

                /* 对角线：物品 id1 单独出现的权重 */
                pairweights[id1][id1] += solval;

                /* 下三角：物品id1 与 id2 (id2 > id1)同箱的LP权重*/
                for (int j = i + 1; j < nconsids; ++j)
                {
                    const int id2 = consids[j];
                    assert(id1 < id2);

                    pairweights[id2][id1] += solval;

                    assert(SCIPisFeasGE(scip, pairweights[id2][id1], 0.0));
                }
            }
        }

        /* 选择 min(p, 1-p) 最大 （最接近0.5）的物品对 */
        SCIP_Real bestvalue = 0.0;
        int id1 = -1;
        int id2 = -1;
        for (int i = 0; i < nitems; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                const SCIP_Real value = std::min(pairweights[i][j], 1.0 - pairweights[i][j]);
                if (bestvalue < value)
                {
                    bestvalue = value;
                    id1 = j;
                    id2 = i;
                }
            }
        }

        assert(SCIPisFeasPositive(scip, bestvalue));
        assert(id1 >= 0 && id1 < nitems);
        assert(id2 >= 0 && id2 < nitems);

        /* 创建子节点 */
        SCIP_NODE* childsame = nullptr;
        SCIP_NODE* childdiffer = nullptr;
        SCIP_CALL(SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)));
        SCIP_CALL(SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)));

        /* 创建 SAME / DIFFER 约束，分别挂载到两个节点上 */
        SCIP_CONS* conssame = nullptr;
        SCIP_CONS* consdiffer = nullptr;
        SCIP_CALL(SCIPcreateConsSamediff(scip, &conssame, "same", id1, id2, ConsType::Same, childsame, true));
        SCIP_CALL(SCIPcreateConsSamediff(scip, &consdiffer, "differ", id1, id2, ConsType::Differ, childdiffer, true));

        SCIP_CALL(SCIPaddConsNode(scip, childsame, conssame, nullptr));
        SCIP_CALL(SCIPaddConsNode(scip, childdiffer, consdiffer, nullptr));

        SCIP_CALL(SCIPreleaseCons(scip, &conssame));
        SCIP_CALL(SCIPreleaseCons(scip, &consdiffer));

        *result = SCIP_BRANCHED;

        return SCIP_OKAY;
    }

} // namespace binpacking