#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

set(SUITE_FILES ../../extra/test_main.cpp ../../extra/test_suite.hpp)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES limits.cpp Jamfile)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../extra PREFIX "_extra" FILES ${SUITE_FILES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../../src PREFIX "_extra" FILES ../../src/src.cpp)
add_executable(boost_url_limits limits.cpp Jamfile ${SUITE_FILES} ../../src/src.cpp)
target_include_directories(boost_url_limits PRIVATE ../../include ../../extra ../../..)
target_compile_definitions(boost_url_limits PRIVATE
    BOOST_URL_MAX_SIZE=16
    BOOST_URL_NO_LIB=1
)
if (BOOST_URL_FIND_PACKAGE_BOOST)
    target_link_libraries(boost_url_limits PRIVATE Boost::headers)
else()
    target_link_libraries(boost_url_limits PRIVATE
        Boost::align
        Boost::config
        Boost::core
        Boost::optional
        Boost::type_traits
        Boost::system
        Boost::variant2)
endif()
target_link_libraries(boost_url_limits INTERFACE Boost::url)
add_test(NAME boost_url_limits COMMAND boost_url_limits)
add_dependencies(boost_url_all_tests boost_url_limits)
