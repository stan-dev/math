# Boost.Geometry
# Example CMakeLists.txt building the Boost.Geometry with wxWidget example
#
# Copyright (c) 2021-2024 Barend Gehrels, Amsterdam, the Netherlands.

# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.8...3.20)

project(wx_widgets_world_mapper)

add_executable(${PROJECT_NAME} wxwidgets_world_mapper.cpp)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

# Link the wxWidgets libraries to our executable
# Assuming it can be found by CMake in /usr/local
find_package(wxWidgets 3.3 COMPONENTS core base REQUIRED CONFIG)
target_link_libraries(${PROJECT_NAME} wxWidgets::wxWidgets)

# Link the Boost.Geometry library to our executable
# By default, it is assumed to be relative to this directory.
target_include_directories(${PROJECT_NAME} PRIVATE .. ../../../../..)

# If this does not work, or you build from elsewhere
# First set BOOST_ROOT
# Then use find_package(Boost)
# Then use target_link_libraries(${PROJECT_NAME} Boost::geometry)
