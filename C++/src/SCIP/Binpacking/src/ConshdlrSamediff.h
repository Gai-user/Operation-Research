#pragma once

#include "objscip/objscip.h"

namespace binpacking
{

    /**
     * 约束类型： Same -> 两物品必须在同一个箱子
     * Differ -> 两物品不能在同一个箱子
     */
    enum class ConsType
    {
        Differ = 0,
        Same = 1
    };
} // namespace binpacking

/**
 * @brief Same/Differ 约束的数据
 *
 * SCIP 用 sizeof(SCIP_ConsData) 分配约束数据内存，因此这个结构体不能放进命名空间。
 * 它记录一个 Ryan/Foster 分支决策： 一对物品(item1, item2) 以及 Same 或 Differ。
 * 额外字段支撑"增量传播"：
 *   - npropagatedvars：上次传播时已检查的变量数（游标）
 *   - propagated：本次节点是否已传播过
 *   - npropagations：传播次数（调试/统计用）
 *   - node：该决策所在的分支节点
 */
struct SCIP_ConsData
{
    int itemid1;                 // 物品1
    int itemid2;                 // 物品2
    binpacking::ConsType type;   // 约束类型
    int npropagatedvars;          // 已传播变量
    int npropagations;           // 传播次数
    unsigned int propagated : 1; // 本次是否已传播
    SCIP_NODE *node;             // 约束所在节点

    SCIP_ConsData(int id1, int id2, binpacking::ConsType t, SCIP_NODE *nd)
        : itemid1(id1), itemid2(id2), type(t),
          npropagatedvars(0), npropagations(0), propagated(0), node(nd)
    {
    }
};

namespace binpacking
{

    /**
     * @brief Same/Differ 约束处理器
     *
     * 该约束处理器不建LP行，而是通过传播吧违反分支决策的装箱变量局部固定位 0：
     *     - Same -> 恰好只包含两物品之一 ，将变量固定为0；
     *     - Differ -> 同时包含两物品，将变量固定为0.
     *
     * enforcement / check / lock 回调原 C 代码为 NULL，C++ 中必须 override（纯虚函数），
     * 实现为空操作直接返回 SCIP_OKAY。
     */
    class ConshdlrSamediff : public scip::ObjConshdlr
    {
    public:
        explicit ConshdlrSamediff(SCIP *scip);
        ~ConshdlrSamediff() override = default;

        virtual SCIP_DECL_CONSDELETE(scip_delete) override;
        virtual SCIP_DECL_CONSTRANS(scip_trans) override;
        virtual SCIP_DECL_CONSENFOLP(scip_enfolp) override;
        virtual SCIP_DECL_CONSENFOPS(scip_enfops) override;
        virtual SCIP_DECL_CONSCHECK(scip_check) override;
        virtual SCIP_DECL_CONSLOCK(scip_lock) override;
        virtual SCIP_DECL_CONSPROP(scip_prop) override;
        virtual SCIP_DECL_CONSACTIVE(scip_active) override;
        virtual SCIP_DECL_CONSDEACTIVE(scip_deactive) override;
        virtual SCIP_DECL_CONSPRINT(scip_print) override;
    };

    /** 创建 SAME/DIFFER 约束（原 SCIPcreateConsSamediff） */
    SCIP_RETCODE SCIPcreateConsSamediff(
        SCIP *scip,
        SCIP_CONS **cons,
        const char *name,
        int itemid1,
        int itemid2,
        ConsType type,
        SCIP_NODE *node,
        SCIP_Bool local);

    int SCIPgetItemid1Samediff(SCIP_CONS *cons);
    int SCIPgetItemid2Samediff(SCIP_CONS *cons);
    ConsType SCIPgetTypeSamediff(SCIP_CONS *cons);
} // namespace binpacking