# Boost.Geometry
#
# Copyright (c) 2024-2024 Barend Gehrels, Amsterdam, the Netherlands.

# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.8...3.20)
project(doxygen_xml2qbk LANGUAGES CXX)

# Find Boost. The actual version is not relevant. It should find it somewhere. If not, set BOOST_ROOT
find_package(Boost 1.81.0 REQUIRED COMPONENTS program_options)

add_executable(${PROJECT_NAME}
    doxygen_xml2qbk.cpp
)

target_compile_features(${PROJECT_NAME}
  PUBLIC 
    cxx_std_14)

# From Boost it uses program_options
target_link_libraries(${PROJECT_NAME}
  PRIVATE
    Boost::program_options
)

# RapidXML is included in the source tree
target_include_directories(${PROJECT_NAME}
  PRIVATE
    .
    contrib/rapidxml-1.13
)
