set(SRC_FILE BranchAndBound/TSP.cpp)
add_library(operationsearch STATIC ${SRC_FILE})
target_compile_features(operationsearch PUBLIC cxx_std_17)

add_subdirectory(SCIP)

install(TARGETS operationsearch
        ARCHIVE DESTINATION lib)