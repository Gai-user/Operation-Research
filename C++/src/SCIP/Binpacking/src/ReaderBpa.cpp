#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "ReaderBpa.h"
#include "scip/scip.h"

#include "ProbDataBinpacking.h"

#define READER_NAME "bpareader"
#define READER_DESC "file reader for binpacking data format"
#define READER_EXTENSION "bpa"

namespace binpacking
{

    SCIP_DECL_READERREAD(ReaderBpa::scip_read)
    {
        assert(scip != nullptr);
        assert(filename != nullptr);
        assert(result != nullptr);

        *result = SCIP_DIDNOTRUN;

        /* 用 std::ifstream 打开文件 */
        std::ifstream file(filename);

        if (!file.is_open())
        {
            SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
            SCIPprintSysError(filename);
            return SCIP_NOFILE;
        }

        /* 第1行，读取问题名 */
        std::string name;
        if (!std::getline(file, name))
        {
            SCIPwarningMessage(scip, "invalid input line 1 in file <%s>: empty file\n", filename);
            return SCIP_READERROR;
        }

        /* 第 2 行：容量、物品数、已知可行解值 */
        std::string line;
        if (!std::getline(file, line))
        {
            return SCIP_READERROR;
        }

        int capacity = 0;
        int nitems = 0;
        int bestsolvalue = 0;
        const int nread = std::sscanf(line.c_str(), "%d %d %d", &capacity, &nitems, &bestsolvalue);

        if (nread < 2)
        {
            SCIPwarningMessage(scip, "invalid input line 2 in file <%s>: <%s>\n", filename, line.c_str());
            return SCIP_READERROR;
        }

        /* 后续每行一个重量；ids 取 0..n-1 */
        std::vector<int> ids;
        std::vector<SCIP_Longint> weights;
        ids.reserve(nitems);
        weights.reserve(nitems);

        SCIP_Bool error = false;
        int nweights = 0;
        int lineno = 2;

        while (!file.eof() && !error)
        {
            if (!std::getline(file, line))
            {
                break;
            }
            ++lineno;

            int weight = 0;

            if (std::sscanf(line.c_str(), "%d", &weight) == 0)
            {
                SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, line.c_str());
                error = TRUE;
                break;
            }
            ids.push_back(nweights);
            weights.push_back(weight);
            ++nweights;

            if (nweights == nitems)
            {
                break;
            }
        }

        if (nweights < nitems)
        {
            SCIPwarningMessage(scip, "set nitems from <%d> to <%d> since the file <%s> only contains <%d> weights\n",
                               nitems, nweights, filename, nweights);
            nitems = nweights;
        }

        if (error)
        {
            return SCIP_READERROR;
        }

        /* 创建问题数据并激活定价器 */
        SCIP_CALL(SCIPprobdataCreate(scip, name.c_str(), ids, weights, static_cast<SCIP_Longint>(capacity)));

        *result = SCIP_SUCCESS;

        return SCIP_OKAY;
    }

    ReaderBpa::ReaderBpa(SCIP *scip)
        : ObjReader(scip, READER_NAME, READER_DESC, READER_EXTENSION)
    {
    }

} // namespace binpacking