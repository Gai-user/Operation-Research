#include <iostream>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "BranchRyanFoster.h"
#include "ConshdlrSamediff.h"
#include "ReaderBpa.h"

/** SCIP 实例的 RAII 生命周期管理器 */
class ScipManager
{
public:
    ScipManager()
    {
        if (SCIPcreate(&scip_) != SCIP_OKAY) {
            scip_ = nullptr;
        }
    }

    ~ScipManager()
    {
        if (scip_ != nullptr) {
            SCIPfree(&scip_);
            BMScheckEmptyMemory();
        }
    }

    ScipManager(const ScipManager&) = delete;
    ScipManager& operator=(const ScipManager&) = delete;
    SCIP* get() noexcept { return scip_; }
    bool valid() const noexcept {return scip_ != nullptr; }
private:
    SCIP* scip_ = nullptr;
};

/*
 * @param argc 命令行参数数量
 * @param argv 命令行参数数组
 * @return SCIP_OKAY 成功，否则返回错误码
 */
static SCIP_RETCODE runSCIP(int argc, char **argv)
{
    ScipManager scipManager;
    if(!scipManager.valid()) {
        std::cerr << "Failed to create SCIP instance" << std::endl;
        return SCIP_ERROR;
    }

    SCIP* scip = scipManager.get();

    /* 启用调试解 */
    SCIPenableDebugSol(scip);

    /**
     * 注册自定义插件
     * 文件读取器
     * 分支规则
     * 约束处理器
     * 定价器不在此注册，在SCIPprobdataCreate 读入数据时通过 SCIPpricerBinpackingActivate 创建并激活
     */
    SCIP_CALL(SCIPincludeObjReader(scip, new binpacking::ReaderBpa(scip), true));
    SCIP_CALL(SCIPincludeObjBranchrule(scip, new binpacking::BranchRyanFoster(scip), true));
    SCIP_CALL(SCIPincludeObjConshdlr(scip, new binpacking::ConshdlrSamediff(scip), true));

    /* 注册默认插件，包括SCIP自带的分支规则、lp求解器等 */
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));

    /* 列生成禁用重启 */
     SCIPsetIntParam(scip, "presolving/maxrestarts", 0); 

    /* 关闭所有割平面算法 */
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, true));

    SCIP_CALL(SCIPprocessShellArguments(scip, argc, argv, "scip.set"));

    return SCIP_OKAY;
}

/**
 * @brief 求解器入口
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