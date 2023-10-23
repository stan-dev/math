#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 DMitry Arkhipov (grisumbras@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

cmake_minimum_required(VERSION 3.8...3.16)

set(BOOST_URL_VERSION 2)
if(BOOST_SUPERPROJECT_VERSION)
    set(BOOST_URL_VERSION ${BOOST_SUPERPROJECT_VERSION})
endif()

project(boost_url VERSION "${BOOST_URL_VERSION}" LANGUAGES CXX)


set(BOOST_URL_IS_ROOT OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(BOOST_URL_IS_ROOT ON)
endif()

if(BOOST_URL_IS_ROOT)
    include(CTest)
endif()

if(NOT BOOST_SUPERPROJECT_VERSION)
    option(BOOST_URL_INSTALL "Install boost::url files" ${BOOST_URL_IS_ROOT})
    option(BOOST_URL_BUILD_TESTS "Build boost::url tests" ${BUILD_TESTING})
    option(BOOST_URL_BUILD_EXAMPLES "Build boost::url examples" ${BOOST_URL_IS_ROOT})
else()
    set(BOOST_URL_BUILD_TESTS ${BUILD_TESTING})
endif()

if (BOOST_SUPERPROJECT_VERSION)
    set(BOOST_URL_FIND_PACKAGE_BOOST OFF)
elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../CMakeLists.txt" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../Jamroot" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../boost-build.jam" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../bootstrap.sh" AND
        EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../libs")
    set(BOOST_URL_FIND_PACKAGE_BOOST OFF)
else()
    set(BOOST_URL_FIND_PACKAGE_BOOST ON)
endif()

if (BOOST_URL_FIND_PACKAGE_BOOST)
    find_package(Boost 1.78.0 REQUIRED COMPONENTS container)
elseif (BOOST_URL_IS_ROOT)
    set(BOOST_URL_UNIT_TEST_LIBRARIES container filesystem unordered json regex)
    set(BOOST_INCLUDE_LIBRARIES url json ${BOOST_URL_UNIT_TEST_LIBRARIES})
    set(BOOST_EXCLUDE_LIBRARIES url)
    set(CMAKE_FOLDER Dependencies)
    add_subdirectory(../.. Dependencies/boost EXCLUDE_FROM_ALL)
    unset(CMAKE_FOLDER)
endif()

include(GNUInstallDirs)
function(boost_url_setup_properties target)
    target_compile_features(${target} PUBLIC cxx_constexpr)
    target_compile_definitions(${target} PUBLIC BOOST_URL_NO_LIB=1)

    if(BOOST_SUPERPROJECT_VERSION)
        target_include_directories(${target} PUBLIC "${PROJECT_SOURCE_DIR}/include")
    else()
        target_include_directories(${target}
            PUBLIC
                "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
                "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )
    endif()

    if (BOOST_URL_FIND_PACKAGE_BOOST)
        target_link_libraries(${target}
                PUBLIC
                Boost::headers
                )
    else()
        target_link_libraries(${target}
            PUBLIC
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
    endif()
endfunction()


file(GLOB_RECURSE BOOST_URL_HEADERS CONFIGURE_DEPENDS
    include/boost/*.hpp
    include/boost/*.ipp
    include/boost/*.natvis
)

set(BOOST_URL_SOURCES src/src.cpp)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/include/boost PREFIX "" FILES ${BOOST_URL_HEADERS})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src PREFIX "" FILES ${BOOST_URL_SOURCES})


add_library(boost_url ${BOOST_URL_HEADERS} ${BOOST_URL_SOURCES})
add_library(Boost::url ALIAS boost_url)
boost_url_setup_properties(boost_url)

if(BUILD_SHARED_LIBS)
    target_compile_definitions(boost_url PUBLIC BOOST_URL_DYN_LINK=1)
else()
    target_compile_definitions(boost_url PUBLIC BOOST_URL_STATIC_LINK=1)
endif()


if(BOOST_URL_INSTALL AND NOT BOOST_SUPERPROJECT_VERSION)
    install(TARGETS boost_url
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
endif()


if(BOOST_URL_BUILD_TESTS)
    add_subdirectory(test)
endif()

if(BOOST_URL_BUILD_EXAMPLES)
    add_subdirectory(example)
endif()
