# Copyright 2018, 2019 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.20)

project(cmake_subdir_test LANGUAGES CXX)

add_subdirectory(../.. boostorg/function_types)

set(deps

# Primary dependencies

config
core
detail
mpl
preprocessor
type_traits

# Secondary dependencies

assert
static_assert
throw_exception
predef
utility
io
)

foreach(dep IN LISTS deps)

  add_subdirectory(../../../${dep} boostorg/${dep})

endforeach()

# --target check

add_executable(main main.cpp)
target_link_libraries(main Boost::function_types)

enable_testing()
add_test(main main)

add_custom_target(check VERBATIM COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
