cmake_minimum_required(VERSION 3.22)

project(BinPacking)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

add_executable(BinPacking
	src/main.cpp
	src/ReaderBpa.cpp
	src/VarDataBinpacking.cpp
	src/ProbDataBinpacking.cpp
	src/PricerBinpacking.cpp
	src/ConshdlrSamediff.cpp
	src/BranchRyanFoster.cpp
)

target_include_directories(BinPacking PRIVATE /workspace/Operation-Research/C++/include/debug/include)
target_link_directories(BinPacking PRIVATE /workspace/Operation-Research/C++/include/debug/lib)

target_link_libraries(BinPacking PRIVATE scip)

