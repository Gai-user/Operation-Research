set(SRC_FILE BranchAndBound/TSP.cpp)

add_library(operationsearch STATIC ${SRC_FILE})

install(TARGETS operationsearch
        ARCHIVE DESTINATION lib)