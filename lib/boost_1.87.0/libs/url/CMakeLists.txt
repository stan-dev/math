#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 DMitry Arkhipov (grisumbras@gmail.com)
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

#-------------------------------------------------
#
# Project
#
#-------------------------------------------------
cmake_minimum_required(VERSION 3.8)
set(BOOST_URL_VERSION 2)
if (BOOST_SUPERPROJECT_VERSION)
    set(BOOST_URL_VERSION ${BOOST_SUPERPROJECT_VERSION})
endif ()
project(boost_url VERSION "${BOOST_URL_VERSION}" LANGUAGES CXX)
set(BOOST_URL_IS_ROOT OFF)
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(BOOST_URL_IS_ROOT ON)
endif ()
set(__ignore__ ${CMAKE_C_COMPILER})

#-------------------------------------------------
#
# Options
#
#-------------------------------------------------
option(BOOST_URL_BUILD_TESTS "Build boost::url tests even if BUILD_TESTING is OFF" OFF)
option(BOOST_URL_BUILD_FUZZERS "Build boost::url fuzzers" OFF)
option(BOOST_URL_BUILD_EXAMPLES "Build boost::url examples" ${BOOST_URL_IS_ROOT})
option(BOOST_URL_MRDOCS_BUILD "Build the target for MrDocs: see mrdocs.yml" OFF)
option(BOOST_URL_DISABLE_THREADS "Disable threads" OFF)
option(BOOST_URL_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

# Check if environment variable BOOST_SRC_DIR is set
if (NOT DEFINED BOOST_SRC_DIR AND DEFINED ENV{BOOST_SRC_DIR})
    set(DEFAULT_BOOST_SRC_DIR "$ENV{BOOST_SRC_DIR}")
else ()
    set(DEFAULT_BOOST_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
endif ()
set(BOOST_SRC_DIR ${DEFAULT_BOOST_SRC_DIR} CACHE STRING "Boost source dir to use when running CMake from this directory")

#-------------------------------------------------
#
# Boost modules
#
#-------------------------------------------------
# The boost super-project requires one explicit dependency per-line.
set(BOOST_URL_DEPENDENCIES
        Boost::align
        Boost::assert
        Boost::config
        Boost::core
        Boost::mp11
        Boost::optional
        Boost::static_assert
        Boost::system
        Boost::throw_exception
        Boost::type_traits
        Boost::variant2
        )
foreach (BOOST_URL_DEPENDENCY ${BOOST_URL_DEPENDENCIES})
    if (BOOST_URL_DEPENDENCY MATCHES "^[ ]*Boost::([A-Za-z0-9_]+)[ ]*$")
        list(APPEND BOOST_URL_INCLUDE_LIBRARIES ${CMAKE_MATCH_1})
    endif ()
endforeach ()
# Conditional dependencies
if (NOT BOOST_URL_MRDOCS_BUILD)
    if (BUILD_TESTING OR BOOST_URL_BUILD_TESTS)
        set(BOOST_URL_UNIT_TEST_LIBRARIES container filesystem unordered)
    endif()
    if (BOOST_URL_BUILD_EXAMPLES)
        set(BOOST_URL_EXAMPLE_LIBRARIES json regex beast)
    endif()
endif()
# Complete dependency list
set(BOOST_INCLUDE_LIBRARIES ${BOOST_URL_INCLUDE_LIBRARIES} ${BOOST_URL_UNIT_TEST_LIBRARIES} ${BOOST_URL_EXAMPLE_LIBRARIES})
set(BOOST_EXCLUDE_LIBRARIES url)

#-------------------------------------------------
#
# Add Boost Subdirectory
#
#-------------------------------------------------
if (BOOST_URL_IS_ROOT)
    set(CMAKE_FOLDER Dependencies)
    # Find absolute BOOST_SRC_DIR
    if (NOT IS_ABSOLUTE ${BOOST_SRC_DIR})
        set(BOOST_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${BOOST_SRC_DIR}")
    endif ()

    # Validate BOOST_SRC_DIR
    set(BOOST_SRC_DIR_IS_VALID ON)
    foreach (F "CMakeLists.txt" "Jamroot" "boost-build.jam" "bootstrap.sh" "libs")
        if (NOT EXISTS "${BOOST_SRC_DIR}/${F}")
            message(STATUS "${BOOST_SRC_DIR}/${F} does not exist. Fallback to find_package.")
            set(BOOST_SRC_DIR_IS_VALID OFF)
            break()
        endif ()
    endforeach ()

    # Create Boost interface targets
    if (BOOST_SRC_DIR_IS_VALID)
        # From BOOST_SRC_DIR
        if (BUILD_SHARED_LIBS)
            set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
        endif()
        set(BOOST_EXCLUDE_LIBRARIES ${PROJECT_NAME})
        set(PREV_BUILD_TESTING ${BUILD_TESTING})
        set(BUILD_TESTING OFF CACHE BOOL "Build the tests." FORCE)
        add_subdirectory(${BOOST_SRC_DIR} Dependencies/boost EXCLUDE_FROM_ALL)
        set(BUILD_TESTING ${PREV_BUILD_TESTING} CACHE BOOL "Build the tests." FORCE)
        set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${BOOST_SRC_DIR}/tools/cmake/include")
    else ()
        # From Boost Package
        find_package(Boost REQUIRED COMPONENTS container)
        foreach (BOOST_INCLUDE_LIBRARY ${BOOST_INCLUDE_LIBRARIES})
            if (NOT TARGET Boost::${BOOST_INCLUDE_LIBRARY})
                add_library(Boost::${BOOST_INCLUDE_LIBRARY} ALIAS Boost::headers)
            endif ()
        endforeach ()
    endif ()
    unset(CMAKE_FOLDER)
endif ()

#-------------------------------------------------
#
# Library
#
#-------------------------------------------------
file(GLOB_RECURSE BOOST_URL_HEADERS CONFIGURE_DEPENDS include/boost/*.hpp include/boost/*.natvis)
file(GLOB_RECURSE BOOST_URL_SOURCES CONFIGURE_DEPENDS src/*.cpp)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/include/boost PREFIX "" FILES ${BOOST_URL_HEADERS})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src PREFIX "url" FILES ${BOOST_URL_SOURCES})
function(boost_url_setup_properties target)
    target_compile_features(${target} PUBLIC cxx_constexpr)
    target_compile_definitions(${target} PUBLIC BOOST_URL_NO_LIB=1)
    if (BOOST_URL_DISABLE_THREADS)
        target_compile_definitions(${target} PUBLIC BOOST_URL_DISABLE_THREADS=1)
    endif()
    target_include_directories(${target} PUBLIC "${PROJECT_SOURCE_DIR}/include")
    target_link_libraries(${target} PUBLIC ${BOOST_URL_DEPENDENCIES})
    target_compile_definitions(${target} PUBLIC $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,BOOST_URL_DYN_LINK=1,BOOST_URL_STATIC_LINK=1>)
    target_compile_definitions(${target} PRIVATE BOOST_URL_SOURCE)
    # issue #868
    # target_compile_options(${target} PUBLIC $<$<CXX_COMPILER_ID:MSVC>:/EHsc>)
endfunction()

if (BOOST_URL_MRDOCS_BUILD)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/mrdocs.cpp" "#include <boost/url.hpp>\n")
    add_library(boost_url_mrdocs "${CMAKE_CURRENT_BINARY_DIR}/mrdocs.cpp")
    boost_url_setup_properties(boost_url_mrdocs)
    target_compile_definitions(boost_url_mrdocs PUBLIC BOOST_URL_MRDOCS BOOST_URL_NO_SSE2)
    return()
endif()

add_library(boost_url ${BOOST_URL_HEADERS} ${BOOST_URL_SOURCES})
add_library(Boost::url ALIAS boost_url)
boost_url_setup_properties(boost_url)


#-------------------------------------------------
#
# Tests
#
#-------------------------------------------------
if (BUILD_TESTING OR BOOST_URL_BUILD_TESTS)
    if (BOOST_URL_IS_ROOT)
        include(CTest)
    endif ()
    add_subdirectory(test)
endif ()

#-------------------------------------------------
#
# Examples
#
#-------------------------------------------------
if (BOOST_URL_BUILD_EXAMPLES)
    add_subdirectory(example)
endif ()


