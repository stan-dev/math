# Copyright 2018 Mike Dev
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt
#
# NOTE: CMake support for Boost.Tuple is currently experimental at best
#       and the interface is likely to change in the future

cmake_minimum_required(VERSION 3.5)
project(BoostTuple)

add_library(boost_tuple INTERFACE)
add_library(Boost::tuple ALIAS boost_tuple)

target_include_directories(boost_tuple INTERFACE include)

target_link_libraries(boost_tuple
    INTERFACE
        Boost::config
        Boost::core
        Boost::static_assert
        Boost::type_traits
)
