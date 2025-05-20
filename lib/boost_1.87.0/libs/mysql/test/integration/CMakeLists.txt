#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

add_executable(
    boost_mysql_integrationtests

    # Utilities
    src/server_features.cpp
    src/metadata_validator.cpp
    src/spotchecks_helpers.cpp
    src/utils.cpp

    # Actual tests
    test/spotchecks.cpp
    test/execution_requests.cpp
    test/crud.cpp
    test/handshake.cpp
    test/prepared_statements.cpp
    test/stored_procedures.cpp
    test/multi_queries.cpp
    test/static_interface.cpp
    test/reconnect.cpp
    test/any_connection.cpp
    test/character_set_tracking.cpp
    test/pipeline.cpp
    test/connection_pool.cpp
    test/db_specific.cpp
    test/database_types.cpp

    # Snippets
    test/snippets/tutorials.cpp
    test/snippets/overview.cpp
    test/snippets/connection_establishment.cpp
    test/snippets/dynamic.cpp
    test/snippets/static.cpp
    test/snippets/prepared_statements.cpp
    test/snippets/multi_resultset.cpp
    test/snippets/multi_function.cpp
    test/snippets/metadata.cpp
    test/snippets/charsets.cpp
    test/snippets/time_types.cpp
    test/snippets/connection_pool.cpp
    test/snippets/sql_formatting.cpp
    test/snippets/sql_formatting_custom.cpp
    test/snippets/pipeline.cpp
    test/snippets/templated_connection.cpp
)
target_include_directories(
    boost_mysql_integrationtests
    PRIVATE
    include
)
target_link_libraries(
    boost_mysql_integrationtests
    PRIVATE
    boost_mysql_testing
)
boost_mysql_common_target_settings(boost_mysql_integrationtests)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.16)
    target_precompile_headers(
        boost_mysql_integrationtests
        PRIVATE
        pch.hpp
    )
endif()

add_test(
    NAME boost_mysql_integrationtests
    COMMAND boost_mysql_integrationtests
)
