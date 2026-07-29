cmake_minimum_required(VERSION 3.22)
project(TSP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(sciptsp
   src/GomoryHuTree.cpp
   src/ProbDataTSP.cpp
   src/ConshdlrSubtour.cpp
   src/ReaderTSP.cpp
   src/main.cpp
)
set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

target_include_directories(sciptsp PRIVATE "/workspace/Operation-Research/SCIPLearn/include/debug/include")
target_link_directories(sciptsp PRIVATE "/workspace/Operation-Research/SCIPLearn/include/debug/lib")
target_link_libraries(sciptsp PRIVATE scip)