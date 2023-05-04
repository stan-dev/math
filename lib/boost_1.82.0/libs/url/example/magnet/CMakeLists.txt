#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

source_group("" FILES
        magnet.cpp
        filter_view.hpp
        )

add_executable(magnet
        magnet.cpp
        filter_view.hpp
        )

set_property(TARGET magnet PROPERTY FOLDER "Examples")
target_link_libraries(magnet PRIVATE Boost::url)
