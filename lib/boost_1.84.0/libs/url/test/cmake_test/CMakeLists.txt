#
# Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#

cmake_minimum_required(VERSION 3.5...3.16)

project(cmake_subdir_test LANGUAGES CXX)
set(__ignore__ ${CMAKE_C_COMPILER})
set(__ignore__ ${CMAKE_C_FLAGS})

if(BOOST_CI_INSTALL_TEST)
    find_package(Boost CONFIG REQUIRED COMPONENTS url)
else()
    set(BOOST_URL_BUILD_TESTS OFF CACHE BOOL "Build the tests." FORCE)
    add_subdirectory(../.. boostorg/url)
endif()

add_executable(main main.cpp)
target_link_libraries(main Boost::url)

if (BUILD_TESTING)
    enable_testing()
    add_test(NAME main COMMAND main)
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
endif()
