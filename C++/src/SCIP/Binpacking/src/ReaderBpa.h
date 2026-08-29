#pragma once

#include "objscip/objscip.h"

namespace binpacking {

/**
 * @brief bpa 文件读取器（对象继承 scip::ObjReader）
 *
 * 解析装箱数据文件（扩展名 .bpa），格式：
 *   - 第 1 行：实例名
 *   - 第 2 行：容量 capacity、物品数 nitems、已知可行解值（三者以空格分隔）
 *   - 之后每行：一个物品的重量
 *
 * 解析完成后调用 SCIPprobdataCreate 把数据注册到 SCIP 并激活定价器。
 * 与原 C 代码 reader_bpa.c 的唯一差别是文件 IO 使用 std::ifstream 而非 SCIPfopen/SCIPfgets。
 */
class ReaderBpa : public scip::ObjReader
{
public:
    ReaderBpa(SCIP* scip);

    ~ReaderBpa() override = default;

    virtual SCIP_DECL_READERREAD(scip_read) override;
};

} // namespace binpacking