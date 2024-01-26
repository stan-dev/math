# Copyright (c) 2023 Klemens D. Morgenstern
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


if(NOT cxx_std_20 IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(STATUS "Boost.Cobalt: not building, compiler doesn't support C++20.")
    return()
endif()

if(MSVC_VERSION AND MSVC_VERSION LESS 1930)
    message(STATUS "Boost.Cobalt: not building, the lowest supported MSVC version is 1930.  ${MSVC_VERSION} is not supported")
    return()
endif()

set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

try_compile(
        BOOST_COBALT_HAS_COROUTINE_INCLUDE
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/coroutine.cpp
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED 20
        OUTPUT_VARIABLE TRY_COMPILE_OUTPUT)

if (NOT BOOST_COBALT_HAS_COROUTINE_INCLUDE)
    message(STATUS "Boost.Cobalt: not building, can't include <coroutine>.")
    message(DEBUG ${TRY_COMPILE_OUTPUT})
    return()
endif()

try_compile(
        BOOST_COBALT_HAS_CONCEPTS
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/concepts.cpp
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED 20
        OUTPUT_VARIABLE TRY_COMPILE_OUTPUT)

if (NOT BOOST_COBALT_HAS_CONCEPTS)
    message(STATUS "Boost.Cobalt: not building, can't include <concepts> or use them.")
    message(DEBUG ${TRY_COMPILE_OUTPUT})
    return()
endif()

set(BOOST_COBALT_SHOULD_USE_CONTAINER OFF)

try_compile(BOOST_COBALT_HAS_STD_PMR
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/memory_resource.cpp
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED 17)

if (NOT BOOST_COBALT_HAS_STD_PMR)
    set(BOOST_COBALT_SHOULD_USE_CONTAINER ON)
endif()

set(BOOST_COBALT_REQUIREMENTS_MATCHED ON)

