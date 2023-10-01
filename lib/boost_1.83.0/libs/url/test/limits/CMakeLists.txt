#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

# Boost.URL library variant for limits
add_library(boost_url_small_limits ${BOOST_URL_HEADERS} ${BOOST_URL_SOURCES})
boost_url_setup_properties(boost_url_small_limits)
target_compile_definitions(boost_url_small_limits PUBLIC BOOST_URL_MAX_SIZE=16 BOOST_URL_NO_LIB=1)

# Test target
add_executable(boost_url_limits limits.cpp Jamfile ${SUITE_FILES})
target_include_directories(boost_url_limits PRIVATE ../../include ../../extra ../../..)
target_link_libraries(boost_url_limits PRIVATE boost_url_small_limits)
set_source_files_properties(limits.cpp PROPERTIES COMPILE_FLAGS ${BOOST_URL_TEST_FLAGS})

# Folders
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES limits.cpp Jamfile)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../extra PREFIX "_extra" FILES ${SUITE_FILES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../src PREFIX "url" FILES ${BOOST_URL_SOURCES})

# CTest target
add_test(NAME boost_url_limits COMMAND boost_url_limits)
add_dependencies(boost_url_all_tests boost_url_limits)
