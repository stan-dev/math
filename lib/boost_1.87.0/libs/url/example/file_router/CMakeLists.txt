#
# Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

add_executable(file_router file_router.cpp)
target_link_libraries(file_router PRIVATE Boost::url Boost::filesystem)
if (TARGET boost_filesystem AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_target_property(FS_IS_IMPORTED boost_filesystem IMPORTED)
    target_compile_options(boost_filesystem PUBLIC $<$<BOOL:${FS_IS_IMPORTED}>:-Wno-error=restrict>)
endif()

source_group("" FILES file_router.cpp)
set_property(TARGET file_router PROPERTY FOLDER "Examples")
