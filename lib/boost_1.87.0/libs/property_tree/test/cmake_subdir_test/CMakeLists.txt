# Copyright 2018-2021 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.20)

project(cmake_subdir_test LANGUAGES CXX)

set(BOOST_INCLUDE_LIBRARIES property_tree)
add_subdirectory(../../../.. boostorg/boost)

add_executable(main main.cpp)
target_link_libraries(main Boost::property_tree)

enable_testing()
add_test(main main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --no-tests=error -C $<CONFIG>)
