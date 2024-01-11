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
    src/get_endpoint.cpp
    src/metadata_validator.cpp
    src/er_network_variant.cpp
    src/sync_errc.cpp
    src/sync_exc.cpp
    src/async_callback.cpp
    src/async_coroutines.cpp
    src/async_coroutinescpp20.cpp

    # Actual tests
    test/spotchecks.cpp
    test/crud.cpp
    test/handshake.cpp
    test/prepared_statements.cpp
    test/stored_procedures.cpp
    test/multi_queries.cpp
    test/static_interface.cpp
    test/reconnect.cpp
    test/db_specific.cpp
    test/database_types.cpp
)
target_include_directories(
    boost_mysql_integrationtests
    PRIVATE
    "include"
)
target_link_libraries(
    boost_mysql_integrationtests
    PRIVATE
    boost_mysql_testing
    Boost::mysql
    Boost::context
)
boost_mysql_common_target_settings(boost_mysql_integrationtests)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.16)
    target_precompile_headers(
        boost_mysql_integrationtests
        PRIVATE
        pch.hpp
    )
endif()

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
