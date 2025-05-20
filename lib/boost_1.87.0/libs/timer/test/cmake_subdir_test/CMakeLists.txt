# Copyright 2018-2021 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.20)

project(cmake_subdir_test LANGUAGES CXX)

# Put boost_timer.dll in the same directory as main.exe
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(../.. boostorg/timer)

set(deps

  # Primary dependencies
  config
  io
  predef
)

foreach(dep IN LISTS deps)

  add_subdirectory(../../../${dep} boostorg/${dep})

endforeach()

add_executable(main main.cpp)
target_link_libraries(main Boost::timer)

enable_testing()
add_test(main main)
