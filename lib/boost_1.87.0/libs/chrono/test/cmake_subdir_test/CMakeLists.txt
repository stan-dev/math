# Copyright 2018, 2019, 2023 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.16)

project(cmake_subdir_test LANGUAGES CXX)

add_subdirectory(../.. boostorg/chrono)

set(deps
# Primary dependencies

assert
config
core
integer
move
mpl
predef
ratio
static_assert
system
throw_exception
type_traits
typeof
utility
winapi

# Secondary dependencies

preprocessor
variant2
io
mp11
)

foreach(dep IN LISTS deps)

  add_subdirectory(../../../${dep} boostorg/${dep})

endforeach()

add_executable(quick ../quick.cpp)
target_link_libraries(quick Boost::chrono)

enable_testing()
add_test(quick quick)

add_custom_target(check VERBATIM COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
