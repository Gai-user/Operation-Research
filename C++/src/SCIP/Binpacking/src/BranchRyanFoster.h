#pragma once

#include "objscip/objscip.h"

namespace binpacking
{

    /**
     * @brief Ryan/Foster 分支规则（对象继承 scip::ObjBranchrule）
     *
     * 针对集合覆盖（列生成）模型，标准变量分支的"取 0 分支"几乎无用——只禁止了指数多个装箱中的
     * 一个，而"取 1 分支"又立刻固定了一大堆物品，导致搜索树极度不平衡。Ryan/Foster 分支的出发点是
     * 选一对物品，强制它们：
     *   - SAME：必须装进同一个箱子；
     *   - DIFFER：禁止装进同一个箱子。
     * 两种情况都允许"两物品都不出现"的装箱存在。
     *
     * ## 物品对的选取
     * 每个装箱变量（列）的 vardata 记录了它出现在哪些 set covering 约束（即包含哪些物品）。
     * 对当前 LP 的所有分数候选变量 v，把它的 LP 值 solval 累加进一个"物品对"三角矩阵 pairweights：
     *   - 对角线 pairweights[id1][id1] += solval：物品 id1 单独出现的 LP 总权重；
     *   - 下三角 pairweights[id2][id1] += solval：物品 id1 与 id2 一起出现的 LP 总权重。
     * 遍历所有候选列后，pairweights[i][j]（i>j）即"物品 j 与 i 同箱"的 LP 分数总需求。
     * 选 min(pairweights[i][j], 1 - pairweights[i][j]) 最大的对——即最接近 0.5 的对，这样两侧分支最平衡。
     *
     * ## 分支如何实现
     * 选出物品对 (id1, id2) 后，在当前节点下创建两个子节点：
     *   - childsame：添加一条 SAME 约束（ConshdlrSamediff）
     *   - childdiffer：添加一条 DIFFER 约束（ConshdlrSamediff）
     * 这些约束本身不建 LP 行，而是由约束处理器在传播时把违反决策的装箱变量局部固定为 0。
     */
    class BranchRyanFoster : public scip::ObjBranchrule
    {
    public:
        BranchRyanFoster(SCIP* scip);
        ~BranchRyanFoster() override = default;

        virtual SCIP_DECL_BRANCHEXECLP(scip_execlp) override;
    };

} // namespace binpacking