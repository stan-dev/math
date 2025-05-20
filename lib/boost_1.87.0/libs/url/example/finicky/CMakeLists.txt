#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

add_executable(finicky finicky.cpp)
target_link_libraries(finicky PRIVATE Boost::url Boost::json Boost::regex)
source_group("" FILES finicky.cpp)
set_property(TARGET finicky PROPERTY FOLDER "Examples")
