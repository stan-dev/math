#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

add_executable(suffix_list suffix_list.cpp)
target_link_libraries(suffix_list PRIVATE Boost::url)
source_group("" FILES suffix_list.cpp)
set_property(TARGET suffix_list PROPERTY FOLDER "Examples")
