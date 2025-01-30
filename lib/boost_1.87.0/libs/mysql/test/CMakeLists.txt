#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Separately compiled Boost.MySQL and Boost.Asio
add_library(
    boost_mysql_compiled
    STATIC
    common/src/boost_asio.cpp
    common/src/boost_mysql.cpp
)
target_compile_definitions(
    boost_mysql_compiled
    PUBLIC
    BOOST_ASIO_SEPARATE_COMPILATION
    BOOST_MYSQL_SEPARATE_COMPILATION
)
target_link_libraries(boost_mysql_compiled PUBLIC boost_mysql)
boost_mysql_common_target_settings(boost_mysql_compiled)

# boost_mysql_testing contains common definitions, includes and settings.
# Note: old versions of cmake require the sources passed to target_sources to be absolute
add_library(boost_mysql_testing INTERFACE)
target_sources(
    boost_mysql_testing
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/common/src/entry_point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/src/utils.cpp
)
target_link_libraries(
    boost_mysql_testing
    INTERFACE
    boost_mysql_compiled
    Boost::unit_test_framework
    Boost::pfr
)
target_include_directories(
    boost_mysql_testing
    INTERFACE
    common/include
)

# Unit testing
add_subdirectory(unit)

# Integration testing
if(BOOST_MYSQL_INTEGRATION_TESTS)
    add_subdirectory(integration)
endif()
