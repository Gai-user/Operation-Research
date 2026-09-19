set(SRC_FILES BranchAndBound/TSP.cpp Simplex/SimplexPrimal/SimplexPrimal.cpp)
add_library(operationsearch STATIC ${SRC_FILES})
target_compile_features(operationsearch PUBLIC cxx_std_17)

install(TARGETS operationsearch
        ARCHIVE DESTINATION lib)