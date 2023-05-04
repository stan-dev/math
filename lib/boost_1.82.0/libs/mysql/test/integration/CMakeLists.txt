#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# This is only processed when we're the main project and
# BOOST_MYSQL_INTEGRATION_TESTS is on

cmake_minimum_required(VERSION 3.12...3.22) # string JOIN
find_package(Boost ${BOOST_MYSQL_VERSION} REQUIRED COMPONENTS context)

add_executable(
    boost_mysql_integrationtests

    # Utilities
    utils/src/get_endpoint.cpp
    utils/src/metadata_validator.cpp
    utils/src/er_network_variant.cpp
    utils/src/sync_errc.cpp
    utils/src/sync_exc.cpp
    utils/src/async_callback.cpp
    utils/src/async_coroutines.cpp
    utils/src/async_coroutinescpp20.cpp

    # Actual tests
    spotchecks.cpp
    crud.cpp
    handshake.cpp
    prepared_statements.cpp
    reconnect.cpp
    db_specific.cpp
    database_types.cpp

    entry_point.cpp
)
target_include_directories(
    boost_mysql_integrationtests
    PRIVATE
    utils/include
)
target_link_libraries(
    boost_mysql_integrationtests
    PRIVATE
    boost_mysql_testing
    Boost::context
)
common_target_settings(boost_mysql_integrationtests)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.16)
    target_precompile_headers(
        boost_mysql_integrationtests
        PRIVATE
        pch.hpp
    )
endif()
set_source_files_properties(
    entry_point.cpp
    PROPERTIES
    SKIP_PRECOMPILE_HEADERS ON
)

# Compose the test filter
if (NOT "$ENV{BOOST_MYSQL_NO_UNIX_SOCKET_TESTS}" STREQUAL "")
    list(APPEND TEST_EXCLUSIONS "!@unix")
endif()
set (DB_SYSTEM $ENV{BOOST_MYSQL_TEST_DB})
if ("${DB_SYSTEM}" STREQUAL "")
    set(DB_SYSTEM "mysql8")
endif()
list(APPEND TEST_EXCLUSIONS "!@skip_${DB_SYSTEM}")
string(JOIN ":" TEST_FILTER ${TEST_EXCLUSIONS})

if ("${TEST_FILTER}" STREQUAL "")
    add_test(
        NAME boost_mysql_integrationtests
        COMMAND boost_mysql_integrationtests
    )
else()
    add_test(
        NAME boost_mysql_integrationtests
        COMMAND boost_mysql_integrationtests "-t" ${TEST_FILTER}
    )
endif()

# If we are using memcheck, then run a subset of the integration tests
# under valgrind. Coroutine tests don't work well under Valgrind, and
# SSL tests are too slow. We do some other exclusions to reduce runtime
if (BOOST_MYSQL_VALGRIND_TESTS)
    string(JOIN ":" TEST_FILTER ${TEST_EXCLUSIONS} "!@ssl" "!@async_coroutine")
    add_memcheck_test(
        NAME boost_mysql_integrationtests_memcheck
        TARGET boost_mysql_integrationtests
        ARGUMENTS "-t" ${TEST_FILTER}
    )
endif()
