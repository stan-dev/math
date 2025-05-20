# Copyright 2022 Alexander Grund
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5...3.16)

project(cmake_subdir_test LANGUAGES CXX)

# Those 2 should work the same
if(DEFINED BOOST_CI_INSTALL_TEST AND BOOST_CI_INSTALL_TEST)
    message("Using installed Boost")
    find_package(boost_locale REQUIRED)
else()
    message("Building Boost")
    set(BOOST_INCLUDE_LIBRARIES locale)
    add_subdirectory(../../../.. deps/boost EXCLUDE_FROM_ALL)
endif()

add_executable(main main.cpp)
target_link_libraries(main Boost::locale Boost::system)

enable_testing()
add_test(NAME main COMMAND main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
