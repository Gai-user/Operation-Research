cmake_minimum_required(VERSION 3.22)
project(VRP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

add_executable(vrp 
               src/main.cpp
               src/PricerVRP.cpp)

target_include_directories(vrp PRIVATE "/workspace/Operation-Research/SCIPLearn/include/debug/include")
target_link_directories(vrp PRIVATE "/workspace/Operation-Research/SCIPLearn/include/debug/lib")

target_link_libraries(vrp PRIVATE scip)