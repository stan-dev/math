#
# Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

add_executable(sanitize sanitize.cpp)
target_link_libraries(sanitize PRIVATE Boost::url)
source_group("" FILES sanitize.cpp)
set_property(TARGET sanitize PROPERTY FOLDER "Examples")
