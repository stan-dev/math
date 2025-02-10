#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 DMitry Arkhipov (grisumbras@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

# Custom target used by the boost super-project
if(NOT TARGET tests)
    add_custom_target(tests)
    set_property(TARGET tests PROPERTY FOLDER Dependencies)
endif()

# Replicate error flags from Jamfile
if (BOOST_URL_WARNINGS_AS_ERRORS)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7)
            set(BOOST_URL_TEST_FLAGS "-Wall -Werror -Wno-unused-but-set-variable -Wno-maybe-uninitialized")
        else()
            set(BOOST_URL_TEST_FLAGS "-Wall -Werror -Wno-unused-but-set-variable")
        endif()
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        set(BOOST_URL_TEST_FLAGS "-Wall -Werror")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT MATCHES "MSVC")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
            set(BOOST_URL_TEST_FLAGS "-Wall -Werror -Wno-unused-but-set-variable")
        else()
            set(BOOST_URL_TEST_FLAGS "-Wall -Werror")
        endif()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" OR CMAKE_CXX_COMPILER_FRONTEND_VARIANT MATCHES "MSVC")
        set(BOOST_URL_TEST_FLAGS "/W4 /WX")
    endif()

    # Print test configuration if running in CI
    # This is useful for debugging CI failures related to warnings which might be false positives
    if (DEFINED ENV{CI})
        message(STATUS "Boost.URL Tests - Compiler ID: ${CMAKE_CXX_COMPILER_ID} / ${CMAKE_CXX_COMPILER_VERSION}")
        if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT)
            message(STATUS "Boost.URL Tests - Compiler Frontend: ${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}")
        endif()
        message(STATUS "Boost.URL Tests - Platform: ${CMAKE_SYSTEM_NAME} / ${CMAKE_SYSTEM_VERSION}")
        message(STATUS "Boost.URL Tests - C++ standard: ${CMAKE_CXX_STANDARD}")
        message(STATUS "Boost.URL Tests - Test error flags: ${BOOST_URL_TEST_FLAGS}")
    endif()
endif()


set(SUITE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/../extra/test_main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../extra/test_suite.hpp)

add_subdirectory(unit)
add_subdirectory(extra)
add_subdirectory(limits)
if (BOOST_URL_BUILD_FUZZERS AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_subdirectory(fuzz)
endif()
