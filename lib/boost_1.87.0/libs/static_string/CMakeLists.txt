#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/static_string
#

cmake_minimum_required(VERSION 3.5...3.16)

project(boost_static_string VERSION "${BOOST_SUPERPROJECT_VERSION}" LANGUAGES CXX)

set(BOOST_STATIC_STRING_IS_ROOT OFF)
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(BOOST_STATIC_STRING_IS_ROOT ON)
endif ()

if (BOOST_STATIC_STRING_IS_ROOT)
    include(CTest)
endif ()

# Options
if (NOT BOOST_SUPERPROJECT_VERSION)
    option(BOOST_STATIC_STRING_INSTALL "Install boost::static_string files" ${BOOST_STATIC_STRING_IS_ROOT})
    option(BOOST_STATIC_STRING_BUILD_TESTS "Build boost::static_string tests" OFF)
else ()
    set(BOOST_STATIC_STRING_BUILD_TESTS ${BUILD_TESTING})
endif ()

# Find boost
if (BOOST_SUPERPROJECT_VERSION)
    set(BOOST_STATIC_STRING_FIND_PACKAGE_BOOST OFF)
elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../CMakeLists.txt" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../Jamroot" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../boost-build.jam" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../bootstrap.sh" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../libs")
    set(BOOST_STATIC_STRING_FIND_PACKAGE_BOOST OFF)
else ()
    set(BOOST_STATIC_STRING_FIND_PACKAGE_BOOST ON)
endif ()

if (BOOST_STATIC_STRING_FIND_PACKAGE_BOOST)
    find_package(Boost 1.78.0 REQUIRED COMPONENTS container)
elseif (BOOST_STATIC_STRING_IS_ROOT)
    set(BOOST_STATIC_STRING_UNIT_TEST_LIBRARIES core)
    set(BOOST_INCLUDE_LIBRARIES static_string assert container_hash static_assert throw_exception utility ${BOOST_STATIC_STRING_UNIT_TEST_LIBRARIES})
    set(BOOST_EXCLUDE_LIBRARIES static_string)
    set(CMAKE_FOLDER Dependencies)
    add_subdirectory(../.. Dependencies/boost EXCLUDE_FROM_ALL)
    unset(CMAKE_FOLDER)
endif ()

# Sources
include(GNUInstallDirs)
file(GLOB_RECURSE BOOST_STATIC_STRING_HEADERS CONFIGURE_DEPENDS
        include/boost/*.hpp
        include/boost/*.ipp
        include/boost/*.natvis
        )

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/include/boost PREFIX "" FILES ${BOOST_STATIC_STRING_HEADERS})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src PREFIX "" FILES ${BOOST_STATIC_STRING_SOURCES})

# Target
add_library(boost_static_string INTERFACE)
add_library(Boost::static_string ALIAS boost_static_string)

target_compile_features(boost_static_string INTERFACE cxx_constexpr)
if (BOOST_SUPERPROJECT_VERSION)
    target_include_directories(boost_static_string INTERFACE "${PROJECT_SOURCE_DIR}/include")
else ()
    target_include_directories(boost_static_string
            INTERFACE
            "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
            )
endif ()

if (BOOST_STATIC_STRING_FIND_PACKAGE_BOOST)
    target_link_libraries(boost_static_string
            INTERFACE
            Boost::headers
            )
else ()
    target_link_libraries(boost_static_string
            INTERFACE
            Boost::assert
            Boost::container_hash
            Boost::core
            Boost::static_assert
            Boost::throw_exception
            Boost::utility
            )
endif ()

if (BOOST_STATIC_STRING_INSTALL AND NOT BOOST_SUPERPROJECT_VERSION)
    install(TARGETS boost_static_string
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            )

    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/boost
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "*.ipp"
            )
endif ()


if (BUILD_TESTING OR BOOST_STATIC_STRING_BUILD_TESTS)
    add_subdirectory(test)
endif ()

