include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://gitee.com/mirrors/googletest.git # 使用gitmee镜像，解决链接不上github的问题
  GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

set(TEST_SOURCES TSPTest.cpp)

add_executable(ortest ${TEST_SOURCES})

target_link_libraries(ortest
                      operationsearch
                      GTest::gtest_main)
gtest_discover_tests(ortest)