# Copyright 2018, 2019, 2021 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.20)

project(cmake_subdir_test LANGUAGES CXX)

add_subdirectory(../.. boostorg/container_hash)

add_subdirectory(../../../config boostorg/config)
add_subdirectory(../../../type_traits boostorg/type_traits)
add_subdirectory(../../../static_assert boostorg/static_assert)
add_subdirectory(../../../describe boostorg/describe)
add_subdirectory(../../../mp11 boostorg/mp11)

add_executable(quick ../quick.cpp)
target_link_libraries(quick Boost::container_hash)

enable_testing()
add_test(quick quick)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
