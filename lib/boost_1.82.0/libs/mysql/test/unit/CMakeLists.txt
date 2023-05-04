#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

add_executable(
    boost_mysql_unittests
    detail/channel/read_buffer.cpp
    detail/channel/message_parser.cpp
    detail/channel/message_reader.cpp
    detail/channel/message_writer_processor.cpp
    detail/channel/message_writer.cpp
    detail/auth/auth_calculator.cpp
    detail/auxiliar/static_string.cpp
    detail/auxiliar/rows_iterator.cpp
    detail/auxiliar/field_type_traits.cpp
    detail/auxiliar/datetime.cpp
    detail/network_algorithms/start_execution_generic.cpp
    detail/network_algorithms/start_query.cpp
    detail/network_algorithms/start_statement_execution.cpp
    detail/network_algorithms/read_some_rows.cpp
    detail/network_algorithms/read_all_rows.cpp
    detail/network_algorithms/query.cpp
    detail/network_algorithms/execute_statement.cpp
    detail/network_algorithms/close_statement.cpp
    detail/network_algorithms/ping.cpp
    detail/protocol/capabilities.cpp
    detail/protocol/null_bitmap_traits.cpp
    detail/protocol/deserialize_errc.cpp
    detail/protocol/serialization_test.cpp
    detail/protocol/deserialize_text_field.cpp
    detail/protocol/deserialize_binary_field.cpp
    detail/protocol/deserialize_row.cpp
    detail/protocol/process_error_packet.cpp
    date.cpp
    datetime.cpp
    field_view.cpp
    field.cpp
    row_view.cpp
    row.cpp
    rows_view.cpp
    rows.cpp
    metadata.cpp
    client_errc.cpp
    common_server_errc.cpp
    mysql_server_errc.cpp
    mariadb_server_errc.cpp
    diagnostics.cpp
    statement.cpp
    execution_state.cpp
    results.cpp
    connection.cpp
    regressions.cpp
    throw_on_error.cpp
    entry_point.cpp
)
target_link_libraries(
    boost_mysql_unittests
    PRIVATE
    boost_mysql_testing
)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.16)
    target_precompile_headers(
        boost_mysql_unittests
        PRIVATE
        pch.hpp
    )
endif()
set_source_files_properties(
    entry_point.cpp
    PROPERTIES
    SKIP_PRECOMPILE_HEADERS ON
)
common_target_settings(boost_mysql_unittests)

if (BOOST_MYSQL_VALGRIND_TESTS)
    add_memcheck_test(
        NAME boost_mysql_unittests_memcheck
        TARGET boost_mysql_unittests
    )
else()
    add_test(
        NAME boost_mysql_unittests
        COMMAND boost_mysql_unittests
    )
    add_dependencies(tests boost_mysql_unittests)
endif()
