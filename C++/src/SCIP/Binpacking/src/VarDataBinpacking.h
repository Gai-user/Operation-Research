#pragma once

#include <vector>

#include "objscip/objscip.h"

/**
 * @brief 装箱问题的变量数据类
 *
 * 每个装箱变量（列）记录其出现在哪些 set covering 约束中（即该装箱包含哪些物品）。
 * SCIP 无法直接从约束 API 反查"变量出现在哪些约束"，所以用 vardata 保存这个"列视图"。
 *
 * 该信息被三处使用：
 *   - BranchRyanFoster：挑选 LP 分数解中"最接近 0.5"的物品对做分支
 *   - ConshdlrSamediff：检查某个装箱是否违反 SAME/DIFFER 决策
 *   - PricerBinpacking：避免重复生成被局部固定为 0 的列
 *
 * 约束 id 数组必须保持升序（`SCIPsortedvecFindInt` 二分查找依赖此不变量）。
 *
 * 通过继承 scip::ObjVardata，用 `SCIPcreateObjVar(..., deleteobject=TRUE)` 创建变量，
 * 变量被 SCIP 释放时由析构函数自动回收 vardata 对象。
 */
namespace binpacking
{

    class VarDataBinpacking : public scip::ObjVardata
    {
    public:
        /** @param consids 约束 ID 数组， 本类内部会排序，因此传入副本即可 */
        explicit VarDataBinpacking(const std::vector<int> &consids);

        /** @return 升序排列的约束 ID 数组 */
        const std::vector<int> &consids() const noexcept { return consids_; }

        /** @return 约束数量 */
        int nconsids() const noexcept { return static_cast<int>(consids_.size()); }

    private:
        std::vector<int> consids_;
    };

    /** 创建装箱变量并绑定 vardata
     *
     * 等价于原 C 代码 SCIPcreateVarBinpacking()：binary 变量、[0,1] 界、按需设置 initial/removable，
     * 并标记变量在预求解中可被删除（SCIPvarMarkDeletable）。
     *
     * @param vardata 可为 nullptr（此时创建无数据的变量）
     */
    SCIP_RETCODE SCIPcreateVarBinpacking(
        SCIP *scip,
        SCIP_VAR **var,
        const char *name,
        SCIP_Real obj,
        SCIP_Bool initial,
        SCIP_Bool removable,
        VarDataBinpacking *vardata);

} // namespace binpacking