#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

add_executable(
    boost_mysql_unittests

    # Helpers
    src/test_stream.cpp
    src/serialization.cpp

    # Actual tests
    test/auth/auth.cpp

    test/protocol/static_buffer.cpp
    test/protocol/capabilities.cpp
    test/protocol/null_bitmap_traits.cpp
    test/protocol/protocol_field_type.cpp
    test/protocol/serialization.cpp
    test/protocol/binary_serialization.cpp
    test/protocol/deserialize_text_field.cpp
    test/protocol/deserialize_binary_field.cpp
    test/protocol/protocol.cpp

    test/channel/read_buffer.cpp
    test/channel/message_parser.cpp
    test/channel/message_reader.cpp
    test/channel/message_writer.cpp
    test/channel/write_message.cpp

    test/execution_processor/execution_processor.cpp
    test/execution_processor/execution_state_impl.cpp
    test/execution_processor/static_execution_state_impl.cpp
    test/execution_processor/results_impl.cpp
    test/execution_processor/static_results_impl.cpp

    test/network_algorithms/read_resultset_head.cpp
    test/network_algorithms/start_execution.cpp
    test/network_algorithms/read_some_rows.cpp
    test/network_algorithms/read_some_rows_dynamic.cpp
    test/network_algorithms/execute.cpp
    test/network_algorithms/close_statement.cpp
    test/network_algorithms/ping.cpp
    test/network_algorithms/read_some_rows_static.cpp

    test/detail/any_stream_impl.cpp
    test/detail/datetime.cpp
    test/detail/row_impl.cpp
    test/detail/rows_iterator.cpp
    test/detail/execution_concepts.cpp
    test/detail/writable_field_traits.cpp
    test/detail/socket_stream.cpp
    test/detail/typing/meta_check_context.cpp
    test/detail/typing/pos_map.cpp
    test/detail/typing/readable_field_traits.cpp
    test/detail/typing/row_traits.cpp
    test/detail/typing/get_type_index.cpp

    test/misc.cpp
    test/multifn.cpp
    test/execution_state.cpp
    test/static_execution_state.cpp
    test/results.cpp
    test/static_results.cpp
    test/resultset_view.cpp
    test/resultset.cpp
    test/client_errc.cpp
    test/common_server_errc.cpp
    test/mysql_server_errc.cpp
    test/mariadb_server_errc.cpp
    test/connection.cpp
    test/date.cpp
    test/datetime.cpp
    test/field_view.cpp
    test/field.cpp
    test/row_view.cpp
    test/row.cpp
    test/rows_view.cpp
    test/rows.cpp
    test/metadata.cpp
    test/diagnostics.cpp
    test/statement.cpp
    test/throw_on_error.cpp
)
target_include_directories(
    boost_mysql_unittests
    PRIVATE
    "include"
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
boost_mysql_common_target_settings(boost_mysql_unittests)

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
