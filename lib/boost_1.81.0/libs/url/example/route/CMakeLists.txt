#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

source_group("" FILES
        route.cpp
        )

add_executable(route
        route.cpp
        )

set_property(TARGET route PROPERTY FOLDER "Examples")
target_link_libraries(route PRIVATE Boost::url Boost::filesystem)
if (TARGET boost_filesystem AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_target_property(FS_IS_IMPORTED boost_filesystem IMPORTED)
    if (FS_IS_IMPORTED)
        target_compile_options(boost_filesystem PUBLIC -Wno-error=restrict)
    endif()
endif()

