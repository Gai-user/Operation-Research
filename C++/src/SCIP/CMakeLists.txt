set(SCIP_ROOT_DIR ${CMAKE_SOURCE_DIR}/include/scip)

if(NOT EXISTS ${SCIP_ROOT_DIR})
    message(FATAL_ERROR "SCIP root directory not found at: ${SCIP_ROOT_DIR}\nPlease check the path in src/SCIP/CMakeLists.txt")
endif()

set(SCIP_INCLUDE_DIR ${SCIP_ROOT_DIR}/include)
set(SCIP_LIB_DIR ${SCIP_ROOT_DIR}/lib)

message(STATUS "Using pre-built SCIP from: ${SCIP_ROOT_DIR}")
message(STATUS "Include dir: ${SCIP_INCLUDE_DIR}")
message(STATUS "Lib dir: ${SCIP_LIB_DIR}")

add_executable(tsp_solver 
            TSP/main.cpp
            TSP/GomoryHuTree.cpp
            TSP/ProData.cpp
            TSP/ReaderTSP.cpp)

# 包含头文件
target_include_directories(tsp_solver PRIVATE
    ${SCIP_INCLUDE_DIR}
)

# 链接scip库
target_link_directories(tsp_solver PRIVATE
    ${SCIP_LIB_DIR}
)

target_link_libraries(tsp_solver scip)