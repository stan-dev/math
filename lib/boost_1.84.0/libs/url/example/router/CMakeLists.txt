#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

add_executable(router router.cpp impl/matches.cpp detail/impl/router.cpp)
target_link_libraries(router PRIVATE Boost::url Boost::beast)

source_group("" FILES router.cpp)
set_property(TARGET router PROPERTY FOLDER "Examples")
