#
# Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#

cmake_minimum_required(VERSION 3.5...3.16)

project(cmake_subdir_test LANGUAGES CXX)

if(BOOST_CI_INSTALL_TEST)
    find_package(boost_url REQUIRED)
else()
    set(DEPENDENCIES
        # boostdep --brief url
        # Primary dependencies
        align assert config core mp11 optional static_assert system throw_exception type_traits variant2
        # Secondary dependencies
        detail move predef utility winapi preprocessor io
    )
    foreach (dep ${DEPENDENCIES})
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../../${dep}/CMakeLists.txt")
            add_subdirectory(../../../${dep} boostorg/${dep})
        endif()
    endforeach()

    add_subdirectory(../.. boostorg/url)
endif()

add_executable(main main.cpp)
target_link_libraries(main Boost::url)

enable_testing()
add_test(NAME main COMMAND main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
