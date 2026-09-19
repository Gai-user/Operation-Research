# 指定 CMake 最低版本
cmake_minimum_required(VERSION 3.20)

# 项目名称和语言
project(PrimalSimplex LANGUAGES CXX)

option(BUILD_TESTING "enable GoogleTest" OFF)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置静态库的输出目录
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# 设置可执行文件的输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(PrimalSimplex SimplexPrimal.cpp
                             main.cpp)