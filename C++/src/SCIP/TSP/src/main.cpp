/**@file   examples/TSP/src/cppmain.cpp
 * @brief  使用 SCIP 作为可调用库的 C++ TSP 示例主文件
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 * @details
 * 本文件是 TSP C++ 示例程序的入口点。它：
 *   1. 创建 SCIP 求解器实例（通过 RAII 包装 ScipContext）
 *   2. 注册 TSP 特定的插件：
 *      - ReaderTSP：解析 .tsp 文件
 *      - ConshdlrSubtour：子环消除约束处理器
 *   3. 加载 SCIP 默认插件（分支定界、割平面、启发式等）
 *   4. 处理命令行参数并执行求解
 *
 * ## 使用方式
 *
 * 编译后通过命令行运行：
 * ~~~
 *   ./sciptsp -f <tsp文件路径> [-l <输出文件>] [其他 SCIP 参数]
 * ~~~
 *
 * 例如：
 * ~~~
 *   ./sciptsp -f ../data/gr17.tsp
 * ~~~
 *
 * ## SCIP_CALL 宏
 *
 * SCIP 使用返回码（SCIP_RETCODE）进行错误处理。大多数 SCIP 函数
 * 返回 SCIP_OKAY（成功）或错误码。SCIP_CALL 宏会检查返回值：
 * 如果不是 SCIP_OKAY，则立即将错误码向上传播。这类似于 Rust 的 `?` 操作符。
 *
 * ## RAII 与 SCIP 生命周期
 *
 * ScipContext 使用 RAII 模式管理 SCIP 的生命周期：
 *   - 构造函数调用 SCIPcreate()
 *   - 析构函数调用 SCIPfree() 和 BMScheckEmptyMemory()（内存泄漏检查）
 *   - 禁止拷贝，确保单一所有权
 */

#include <iostream>
#include <memory>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "ReaderTSP.h"
#include "ConshdlrSubtour.h"

/**
 * @brief SCIP 实例的 RAII 生命周期管理器
 *
 * 在构造时创建 SCIP 实例，在析构时自动释放所有资源。
 * 确保即使发生异常也能正确清理 SCIP 分配的内存。
 *
 * 使用模式：
 * ~~~
 *   ScipContext ctx;
 *   if (!ctx.valid()) { // 处理错误 }
 *   SCIP* scip = ctx.get();
 *   // 使用 scip...
 *   // 析构时自动释放
 * ~~~
 */
class ScipContext
{
public:
    ScipContext()
    {
        if (SCIPcreate(&scip_) != SCIP_OKAY)
            scip_ = nullptr;
    }

    ~ScipContext()
    {
        if (scip_ != nullptr)
        {
            SCIPfree(&scip_);
            // 检查是否有未释放的 SCIP 内存块（调试用）
            BMScheckEmptyMemory();
        }
    }

    ScipContext(const ScipContext &) = delete;
    ScipContext &operator=(const ScipContext &) = delete;

    /** @return SCIP 原始指针 */
    SCIP *get() noexcept { return scip_; }
    /** @return SCIP 实例是否成功创建 */
    bool valid() const noexcept { return scip_ != nullptr; }

private:
    SCIP *scip_ = nullptr;
};

/**
 * @brief 初始化 SCIP 求解器并注册 TSP 插件
 *
 * 执行以下步骤：
 *   1. 创建 SCIP 实例（通过 ScipContext RAII 包装）
 *   2. 启用调试解（debug solution）以辅助开发
 *   3. 注册 TSP 文件读取器（ReaderTSP）：识别 .tsp 文件扩展名
 *   4. 注册子环消除约束处理器（ConshdlrSubtour）：处理子环约束
 *   5. 加载 SCIP 默认插件：分支定界规则、割平面生成器、启发式算法等
 *   6. 处理命令行参数：解析 -f 指定输入文件，-l 指定日志文件等
 *      使用默认设置文件 "sciptsp.set"（如果存在）
 *
 * @param argc 命令行参数数量
 * @param argv 命令行参数数组
 * @return SCIP_OKAY 成功，否则返回错误码
 */
static SCIP_RETCODE runSCIP(int argc, char **argv)
{
    ScipContext ctx;
    if (!ctx.valid())
    {
        std::cerr << "Failed to create SCIP instance" << std::endl;
        return SCIP_NOMEMORY;
    }

    SCIP *scip = ctx.get();

    // 启用调试解：允许在求解过程中记录和检查中间解
    SCIPenableDebugSol(scip);

    // 注册 TSP 自定义插件
    // TRUE 参数表示 SCIP 接管对象所有权，在 SCIPfree 时自动释放
    SCIP_CALL(SCIPincludeObjReader(scip, new tsp::ReaderTSP(scip), true));
    SCIP_CALL(SCIPincludeObjConshdlr(scip, new tsp::ConshdlrSubtour(scip), true));

    // 加载 SCIP 默认插件（分支定界、割平面、启发式、节点选择等）
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));

    // 处理命令行参数
    // sciptsp.set 是默认设置文件，可在其中配置求解参数
    SCIP_CALL(SCIPprocessShellArguments(scip, argc, argv, "sciptsp.set"));

    return SCIP_OKAY;
}

/**
 * @brief TSP 求解器入口点
 *
 * 调用 runSCIP() 初始化并求解，处理可能的错误。
 *
 * @param argc 命令行参数数量
 * @param argv 命令行参数数组
 * @return 0 成功，-1 失败
 */
int main(int argc, char** argv)
{
   SCIP_RETCODE retcode = runSCIP(argc, argv);
   if (retcode != SCIP_OKAY)
   {
      SCIPprintError(retcode);
      return -1;
   }
   return 0;
}