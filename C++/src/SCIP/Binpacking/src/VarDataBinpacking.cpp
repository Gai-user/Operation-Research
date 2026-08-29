#include "VarDataBinpacking.h"

#include <algorithm>

namespace binpacking
{

    VarDataBinpacking::VarDataBinpacking(const std::vector<int> &consids) : consids_(consids)
    {
        std::sort(consids_.begin(), consids_.end());
    }

    SCIP_RETCODE SCIPcreateVarBinpacking(
        SCIP *scip,
        SCIP_VAR **var,
        const char *name,
        SCIP_Real obj,
        SCIP_Bool initial,
        SCIP_Bool removable,
        VarDataBinpacking *vardata)
    {
        /* SCIPcreateObjVar 直接设置 initial/removable，并让 SCIP 在变量释放时删除 vardata 对象 */
        SCIP_CALL(SCIPcreateObjVar(scip, var, name, 0.0, 1.0, obj, SCIP_VARTYPE_BINARY, initial, removable, vardata, true));

        SCIPvarMarkDeletable(*var);

        return SCIP_OKAY;
    }
} // namespace binpacking