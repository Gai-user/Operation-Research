#include <iostream>

/* include SCIP components */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* include TSP specific components */
#include "ReaderTSP.h"

using namespace scip;
using namespace tsp;
using namespace std;

SCIP_RETCODE runSCIP(int argc, char** argv)
{
    SCIP* scip = NULL;

    SCIP_CALL(SCIPcreate(&scip));

    SCIPenableDebugSol(scip);

    SCIP_CALL(SCIPincludeObjReader(scip, new ReaderTSP(scip), TRUE));

    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

    SCIP_CALL(SCIPprocessShellArguments(scip, argc, argv, "sciptsp.set"));

    SCIP_CALL(SCIPfree(&scip));

    BMScheckEmptyMemory();
    return SCIP_OKAY;
}

int main(int argc, char** argv) {
    SCIP_RETCODE retcode;

    retcode = runSCIP(argc, argv);

    if (retcode != SCIP_OKAY) {
        SCIPprintError(retcode);
        return -1;
    }

    return 0;
}