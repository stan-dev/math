# Boost.Geometry
# Example CMakeLists.txt building the Boost.Geometry with Qt example
#
# Copyright (c) 2021-2024 Barend Gehrels, Amsterdam, the Netherlands.

# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.16)
project(qt_world_mapper LANGUAGES CXX)

find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)

qt_standard_project_setup()

qt_add_executable(${PROJECT_NAME}
    qt_world_mapper.cpp
    qt_world_mapper.hpp
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    WIN32_EXECUTABLE TRUE
    MACOSX_BUNDLE TRUE
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ..
    ../../../../..
)
