#pragma once

#include <utility>
#include <vector>

#include "objscip/objscip.h"

namespace binpacking
{

    /**
     * @brief 装箱问题数据类（继承 scip::ObjProbData）
     *
     * 保存装箱实例的全部信息：物品 id/重量、箱子容量、所有已生成的变量（列）、每个物品的
     * set covering 约束。所有插件通过 `SCIPgetObjProbData(scip)` 拿到本对象再调用访问器。
     *
     * ## 数据所有权
     *   - vars_ / conss_：SCIP 变量/约束指针。添加时 SCIPcaptureVar/SCIPcaptureCons，
     *     释放时（scip_delorig / scip_deltrans）逐个 SCIPreleaseVar/SCIPreleaseCons。
     *   - 数组本身用 std::vector 管理，无需手工 realloc。
     *
     * ## 列生成联动
     * 问题数据注册了一个 "addedvar" 事件处理器：每当有变量加入 SCIP（SCIPaddVar 或
     * 定价器 SCIPaddPricedVar 触发 SCIP_EVENTTYPE_VARADDED），回调把该变量追加到 vars_。
     *
     * ## 生命周期
     * 对应四个回调：scip_delorig（释放原始数据）、scip_trans（变换到 presolved 空间）、
     * scip_deltrans（释放变换数据）、scip_initsol/scip_exitsol（求解开始/结束时 catch/drop 事件）。
     */
    class ProbDataBinpacking : public scip::ObjProbData
    {
    public:
        /** 构造：保存 ids/weights/capacity 的副本；conss 为每物品 set covering 约束数组（可选，vars 随后由 createInitialColumns 填充） */
        ProbDataBinpacking(std::vector<int> ids, std::vector<SCIP_Longint> weights, SCIP_Longint capacity, std::vector<SCIP_CONS *> conss = {});

        ~ProbDataBinpacking() override = default;

        /* 访问器 */
        const std::vector<int> &ids() const noexcept { return ids_; }
        const std::vector<SCIP_Longint> &weights() const noexcept { return weights_; }
        int nitems() const noexcept { return nitems_; }
        SCIP_Longint capacity() const noexcept { return capacity_; }
        const std::vector<SCIP_VAR *> &vars() const noexcept { return vars_; }
        const std::vector<SCIP_CONS *> &conss() const noexcept { return conss_; }
        int nvars() const noexcept { return static_cast<int>(vars_.size()); }

        /* 将新变量加入 vars_ 并 capture (由 "addedvar" 事件回调调用)*/
        void addVar(SCIP *scip, SCIP_VAR *var);

        /**
         * 设置每个物品的 set covering 约束数组 - 仅在 SCIPprobdataCreate 中使用一次
         * 问题对象在 SCIPcreateObjProb 时先注册，此时 conss 尚未创建， 随后创建约束并一次性挂入
         */
        void setConss(std::vector<SCIP_CONS *> conss) { conss_ = std::move(conss); }

        // ---- objscip 回调 ----
        SCIP_RETCODE scip_delorig(SCIP *scip) override;
        SCIP_RETCODE scip_trans(SCIP *scip, scip::ObjProbData **objprobdata, SCIP_Bool *deleteobject) override;
        SCIP_RETCODE scip_deltrans(SCIP *scip) override;
        SCIP_RETCODE scip_initsol(SCIP *scip) override;
        SCIP_RETCODE scip_exitsol(SCIP *scip, SCIP_Bool restart) override;

    private:
        std::vector<SCIP_VAR *> vars_;      ///< 所有已生成变量（按生成顺序）
        std::vector<SCIP_CONS *> conss_;    ///< 每个物品一个 set covering 约束
        std::vector<SCIP_Longint> weights_; ///< 物品重量
        std::vector<int> ids_;              ///< 物品 id
        int nitems_;                        ///< 物品数
        SCIP_Longint capacity_;             ///< 箱子容量

        /** 释放 vars_ 与 conss_ 中所有指针（scip_delorig / scip_deltrans 共用） */
        void freeVarsConss(SCIP *scip);
    };

    /** 创建装箱问题数据并注册到 SCIP（由 ReaderBpa 调用）
     *
     * 执行：注册 "addedvar" 事件处理器（若无）→ 为每个物品建 set covering 约束（modifiable=TRUE）
     * → 建初始列 → 设目标函数为最小化 + 整数值 → SCIPcreateObjProb 注册本数据对象 → 激活定价器。
     */
    SCIP_RETCODE SCIPprobdataCreate(
        SCIP *scip,
        const char *probname,
        const std::vector<int> &ids,
        const std::vector<SCIP_Longint> &weights,
        SCIP_Longint capacity);
        
} // namespace binpacking