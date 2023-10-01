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
    set_property(TARGET tests PROPERTY FOLDER _deps)
endif()

# Custom target for all tests
if(NOT TARGET boost_url_all_tests)
    add_custom_target(boost_url_all_tests)
    set_property(TARGET boost_url_all_tests PROPERTY FOLDER Dependencies)
endif()

# Replicate error flags from Jamfile
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION_MAJOR EQUAL 7)
    set(BOOST_URL_TEST_FLAGS "-Wall -Werror -Wno-unused-but-set-variable -Wno-maybe-uninitialized")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(BOOST_URL_TEST_FLAGS "-Wall -Werror -Wno-unused-but-set-variable")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION_MAJOR EQUAL 4)
    set(BOOST_URL_TEST_FLAGS "-Wall -Werror -Wno-unused-but-set-variable")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(BOOST_URL_TEST_FLAGS "-Wall -Werror")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(BOOST_URL_TEST_FLAGS "/W4 /WX")
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
