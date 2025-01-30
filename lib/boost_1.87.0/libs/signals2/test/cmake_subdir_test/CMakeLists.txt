# Copyright 2018, 2019 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.20)

project(cmake_subdir_test LANGUAGES CXX)

add_subdirectory(../.. boostorg/signals2)

# boostdep --brief signals2

set(deps

# Primary dependencies

assert
bind
config
core
function
iterator
move
mpl
optional
parameter
preprocessor
smart_ptr
throw_exception
tuple
type_traits
variant

# Secondary dependencies

static_assert
concept_check
detail
function_types
fusion
utility
predef
mp11
container_hash
integer
type_index
describe
functional
typeof
io
)

foreach(dep IN LISTS deps)

  add_subdirectory(../../../${dep} boostorg/${dep})

endforeach()

# --target check

add_executable(main main.cpp)
target_link_libraries(main Boost::signals2)

enable_testing()
add_test(main main)

add_custom_target(check VERBATIM COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
