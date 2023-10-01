#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

source_group("" FILES
        mailto.cpp
        )

add_executable(mailto
        mailto.cpp
        )

set_property(TARGET mailto PROPERTY FOLDER "Examples")
target_link_libraries(mailto PRIVATE Boost::url)
