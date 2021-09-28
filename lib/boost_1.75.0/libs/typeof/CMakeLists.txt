# Copyright 2018 Mike Dev
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt
#
# NOTE: CMake support for Boost.typeof is currently experimental at best
#       and the interface is likely to change in the future

cmake_minimum_required(VERSION 3.5)
project(BoostTypeof LANGUAGES CXX)

add_library(boost_typeof INTERFACE)
add_library(Boost::typeof ALIAS boost_typeof)

target_include_directories(boost_typeof INTERFACE include)

target_link_libraries(boost_typeof
    INTERFACE
        Boost::config
        Boost::preprocessor
        Boost::type_traits
)

