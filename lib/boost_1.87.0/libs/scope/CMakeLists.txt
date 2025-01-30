# Copyright 2023 Andrey Semashev
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# https://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.5...3.20)

project(boost_scope VERSION "${BOOST_SUPERPROJECT_VERSION}" LANGUAGES CXX)

add_library(boost_scope INTERFACE)
add_library(Boost::scope ALIAS boost_scope)

target_include_directories(boost_scope INTERFACE include)

target_link_libraries(boost_scope
    INTERFACE
        Boost::config
        Boost::core
        Boost::type_traits
)

target_compile_features(boost_scope
    INTERFACE
        cxx_constexpr
        cxx_noexcept
        cxx_rvalue_references
        cxx_decltype
        cxx_alias_templates
        cxx_variadic_templates
        cxx_variadic_macros
        cxx_defaulted_functions
        cxx_deleted_functions
        cxx_default_function_template_args
        cxx_delegating_constructors
)

if (BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")
    add_subdirectory(test)
endif()
