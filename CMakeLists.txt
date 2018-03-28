cmake_minimum_required(VERSION 3.9)
project(math)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(lib/gtest_1.7.0)

include_directories(
  ${gtest_SOURCE_DIR}/include
)

include_directories(
  ${CMAKE_SOURCE_DIR}
)


add_executable(broken_test
        test/unit/math/rev/core/agrad_test.cpp
)
