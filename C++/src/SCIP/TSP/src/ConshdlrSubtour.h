/**@file   ConshdlrSubtour.h
 * @brief  TSP 子环消除约束的 C++ 约束处理器
 * @author Timo Berthold
 *
 * @details
 * 在 TSP（旅行商问题）的 SCIP 求解中，约束处理器负责两项关键任务：
 *
 * ## 1. 解的可行性检查 (Check)
 * 检查当前解是否构成完整的哈密顿环（访问所有节点恰好一次）。
 * 如果不是，则解不可行。
 *
 * ## 2. 割平面分离 (Separation)
 * 当 LP 松弛解存在子环时，通过 Gomory-Hu 树算法生成割平面
 * （子环消除不等式），割断非法的子环解，加强 LP 松弛。
 *
 * ## 子环消除不等式的数学形式
 *
 * 对于图的任意真子集 S ⊂ V（0 < |S| < |V|），以下不等式必须成立：
 * ~~~
 *   Σ_{u∈S, v∉S} x_{uv} ≥ 2
 * ~~~
 *
 * 即：跨越割 (S, V\S) 的边数至少为 2。如果 LP 解中该和小于 2，
 * 说明子集 S 在解中形成了孤立子环，需要添加割平面切断它。
 *
 * ## 与其他组件的关系
 *   - **ReaderTSP**：解析 TSPLIB 文件 → 创建 Graph + 度约束 + ConshdlrSubtour
 *   - **ProbDataTSP**：持有 Graph 的所有权，约束处理器通过 SCIP_ConsData 获取引用
 *   - **GomoryHuTree**：提供 ghc_tree() 函数，生成 Gomory-Hu 割
 */

#pragma once

#include "objscip/objscip.h"
#include "ProbDataTSP.h"

/**
 * @brief 子环消除约束的自定义数据
 *
 * SCIP 的 C API 通过 SCIP_ConsData 管理约束数据。
 * 这里存储对 TSP 图的非拥有引用（Graph*），实际所有权由 ProbDataTSP 持有。
 *
 * @note 必须定义在全局命名空间 ::SCIP_ConsData，
 *       因为 SCIP 使用 sizeof(SCIP_ConsData) 来分配约束数据内存块。
 */
struct SCIP_ConsData
{
    Graph *graph; ///< 非拥有引用，指向 TSP 的完全图

    explicit SCIP_ConsData(Graph *g) : graph(g) {}
};

namespace tsp
{

    /**
     * @brief TSP 子环消除约束的 C++ 约束处理器
     *
     * 继承自 scip::ObjConshdlr，实现 SCIP 约束处理器的所有必要回调。
     *
     * ## 回调函数说明
     *
     * | 回调 | 触发时机 | 本处理器的作用 |
     * |------|---------|--------------|
     * | scip_delete | 约束被删除时 | 无操作（数据由 ProbDataTSP 管理） |
     * | scip_trans | 问题变换到 presolved 空间时 | 深拷贝图数据 |
     * | scip_sepalp | LP 松弛求解后 | 生成 Gomory-Hu 割平面分离子环 |
     * | scip_sepasol | 整数可行解找到后 | 检查是否存在子环 |
     * | scip_enfolp | LP 解违反约束且需强制执行时 | 同 sepalp，但强制添加割 |
     * | scip_enfops | 伪解违反约束时 | 检查可行性 |
     * | scip_check | 验证解是否可行时 | 检查是否形成完整哈密顿环 |
     * | scip_prop | 约束传播时 | 在 presolve 后更新图数据 |
     * | scip_lock | 变量锁定 | 为约束涉及的所有变量添加锁 |
     * | scip_delvars | 变量被删除时 | 无操作 |
     * | scip_print | 输出约束到文件时 | 打印图的边信息 |
     *
     * ## 优先级说明
     *
     * 构造函数中设置了较高/较低的优先级值（1000000 / -2000000），
     * 确保分离回调在其他约束处理器之前/之后执行，以优化割平面生成顺序。
     */
    class ConshdlrSubtour : public scip::ObjConshdlr
    {
    public:
        explicit ConshdlrSubtour(SCIP *scip)
            : ObjConshdlr(scip, "subtour", "TSP subtour elimination constraints",
                          1000000, -2000000, -2000000, 1, -1, 1, 0,
                          FALSE, FALSE, TRUE, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_FAST)
        {
        }

        virtual ~ConshdlrSubtour() = default;

        /** @brief 约束被删除时的清理回调 */
        virtual SCIP_DECL_CONSDELETE(scip_delete);
        /** @brief 问题变换时的数据转换回调 */
        virtual SCIP_DECL_CONSTRANS(scip_trans);
        /** @brief LP 分离回调：在 LP 松弛解中寻找违反的割 */
        virtual SCIP_DECL_CONSSEPALP(scip_sepalp);
        /** @brief 整数解分离回调：在整数解中寻找子环 */
        virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);
        /** @brief LP 强制执行回调：强制切断违反 LP 松弛的子环 */
        virtual SCIP_DECL_CONSENFOLP(scip_enfolp);
        /** @brief 伪解强制执行回调：检查伪解的可行性 */
        virtual SCIP_DECL_CONSENFOPS(scip_enfops);
        /** @brief 可行性检查回调：验证解是否构成完整哈密顿环 */
        virtual SCIP_DECL_CONSCHECK(scip_check);
        /** @brief 约束传播回调 */
        virtual SCIP_DECL_CONSPROP(scip_prop);
        /** @brief 变量锁定回调：锁定约束涉及的所有变量 */
        virtual SCIP_DECL_CONSLOCK(scip_lock);
        /** @brief 变量删除回调 */
        virtual SCIP_DECL_CONSDELVARS(scip_delvars);
        /** @brief 约束打印回调：输出约束信息到文件 */
        virtual SCIP_DECL_CONSPRINT(scip_print);

        /** @brief 约束处理器可被克隆 */
        virtual SCIP_DECL_CONSHDLRISCLONEABLE(iscloneable) { return TRUE; }
        /** @brief 克隆回调：深拷贝约束处理器及其数据 */
        virtual SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable *clone);
        /** @brief 约束拷贝回调：在子 SCIP 实例中复制约束 */
        virtual SCIP_DECL_CONSCOPY(scip_copy);
    };

    /**
 * @brief 创建并注册一个 TSP 子环消除约束
 *
 * @param scip       SCIP 实例
 * @param cons       输出：新创建的约束
 * @param name       约束名称
 * @param graph      TSP 图数据（非拥有引用）
 * @param initial    是否为初始约束
 * @param separate   是否启用分离
 * @param enforce    是否启用强制执行
 * @param check      是否启用可行性检查
 * @param propagate  是否启用传播
 * @param local      是否为局部约束
 * @param modifiable 是否可被修改
 * @param dynamic    是否为动态约束
 * @param removable  是否可被移除
 * @return SCIP_OKAY 成功，否则错误码
 */
SCIP_RETCODE SCIPcreateConsSubtour(
   SCIP*                 scip,
   SCIP_CONS**           cons,
   const char*           name,
   Graph*                graph,
   SCIP_Bool             initial,
   SCIP_Bool             separate,
   SCIP_Bool             enforce,
   SCIP_Bool             check,
   SCIP_Bool             propagate,
   SCIP_Bool             local,
   SCIP_Bool             modifiable,
   SCIP_Bool             dynamic,
   SCIP_Bool             removable
   );

} // end tsp
