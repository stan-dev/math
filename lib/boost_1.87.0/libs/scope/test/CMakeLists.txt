# Copyright 2023 Peter Dimov
# Copyright 2023-2024 Andrey Semashev
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

include(BoostTest OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if (NOT HAVE_BOOST_TEST)
    return()
endif()

set(BOOST_TEST_LINK_LIBRARIES Boost::scope)
include_directories(common)

set(BOOST_TEST_COMPILE_FEATURES
    cxx_nullptr
    cxx_lambdas
    cxx_auto_type
    cxx_uniform_initialization
)

if (WIN32)
    set(BOOST_TEST_COMPILE_DEFINITIONS "_CRT_SECURE_NO_WARNINGS" "_CRT_SECURE_NO_DEPRECATE")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(BOOST_TEST_COMPILE_OPTIONS "-Wall" "-Wextra" "-Werror")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(BOOST_TEST_COMPILE_OPTIONS "/W4" "/WX")
endif()

file(GLOB RUN_TESTS LIST_DIRECTORIES OFF CONFIGURE_DEPENDS run/*.cpp)

foreach(TEST IN LISTS RUN_TESTS)
    boost_test(TYPE run SOURCES ${TEST})
endforeach()

unset(BOOST_TEST_COMPILE_OPTIONS)

file(GLOB COMPILE_TESTS LIST_DIRECTORIES OFF CONFIGURE_DEPENDS compile/*.cpp)

foreach(TEST IN LISTS COMPILE_TESTS)
    if("${TEST}" STREQUAL "${CMAKE_CURRENT_LIST_DIR}/compile/self_contained_header.cpp")
        continue()
    endif()

    boost_test(TYPE compile SOURCES ${TEST})
endforeach()

file(GLOB COMPILE_FAIL_TESTS LIST_DIRECTORIES OFF CONFIGURE_DEPENDS compile_fail/*.cpp)

foreach(TEST IN LISTS COMPILE_FAIL_TESTS)
    boost_test(TYPE compile-fail SOURCES ${TEST})
endforeach()
