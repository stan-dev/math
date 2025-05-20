# Copyright 2018 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.16)

project(cmake_subdir_test LANGUAGES CXX)

add_subdirectory(../.. boostorg/mpl)

add_subdirectory(../../../assert boostorg/assert)
add_subdirectory(../../../config boostorg/config)
add_subdirectory(../../../core boostorg/core)
add_subdirectory(../../../io boostorg/io)
add_subdirectory(../../../predef boostorg/predef)
add_subdirectory(../../../preprocessor boostorg/preprocessor)
add_subdirectory(../../../static_assert boostorg/static_assert)
add_subdirectory(../../../throw_exception boostorg/throw_exception)
add_subdirectory(../../../type_traits boostorg/type_traits)
add_subdirectory(../../../utility boostorg/utility)

add_executable(main main.cpp)
target_link_libraries(main Boost::mpl)

enable_testing()
add_test(NAME main COMMAND main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
