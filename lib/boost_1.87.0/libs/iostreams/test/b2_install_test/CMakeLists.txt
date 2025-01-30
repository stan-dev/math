# Copyright 2018, 2019, 2021 Peter Dimov
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.8...3.20)

project(boost_iostreams_install_test LANGUAGES CXX)

if(BOOST_RUNTIME_LINK STREQUAL "static")
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

set(BOOST_HINTS)

if(USE_STAGED_BOOST)
  set(BOOST_HINTS HINTS ../../../../stage)
endif()

find_package(boost_iostreams CONFIG REQUIRED ${BOOST_HINTS})

add_executable(test_gzip test_gzip.cpp)
target_link_libraries(test_gzip Boost::iostreams)
target_compile_features(test_gzip PUBLIC cxx_std_11)

add_executable(test_bzip2 test_bzip2.cpp)
target_link_libraries(test_bzip2 Boost::iostreams)
target_compile_features(test_bzip2 PUBLIC cxx_std_11)

enable_testing()
add_test(NAME test_gzip COMMAND test_gzip WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME test_bzip2 COMMAND test_bzip2 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(check VERBATIM COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --no-tests=error -C $<CONFIG>)
