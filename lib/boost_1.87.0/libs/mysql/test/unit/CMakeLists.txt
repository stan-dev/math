#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

add_executable(
    boost_mysql_unittests

    # Helpers
    src/utils.cpp

    # Actual tests
    test/auth/auth.cpp

    test/protocol/static_buffer.cpp
    test/protocol/capabilities.cpp
    test/protocol/null_bitmap.cpp
    test/protocol/protocol_field_type.cpp
    test/protocol/frame_header.cpp
    test/protocol/serialization_context.cpp
    test/protocol/deserialization_context.cpp
    test/protocol/protocol_types.cpp
    test/protocol/text_protocol.cpp
    test/protocol/binary_protocol.cpp
    test/protocol/serialization.cpp
    test/protocol/deserialization.cpp

    test/sansio/read_buffer.cpp
    test/sansio/message_reader.cpp
    test/sansio/top_level_algo.cpp

    test/sansio/read_resultset_head.cpp
    test/sansio/start_execution.cpp
    test/sansio/read_some_rows.cpp
    test/sansio/read_some_rows_dynamic.cpp
    test/sansio/execute.cpp
    test/sansio/close_statement.cpp
    test/sansio/set_character_set.cpp
    test/sansio/ping.cpp
    test/sansio/reset_connection.cpp
    test/sansio/prepare_statement.cpp
    test/sansio/run_pipeline.cpp

    test/execution_processor/execution_processor.cpp
    test/execution_processor/execution_state_impl.cpp
    test/execution_processor/static_execution_state_impl.cpp
    test/execution_processor/results_impl.cpp
    test/execution_processor/static_results_impl.cpp

    test/connection_pool/sansio_connection_node.cpp
    test/connection_pool/connection_pool_impl.cpp

    test/detail/datetime.cpp
    test/detail/row_impl.cpp
    test/detail/rows_iterator.cpp
    test/detail/execution_concepts.cpp
    test/detail/writable_field_traits.cpp
    test/detail/socket_stream.cpp
    test/detail/connect_params_helpers.cpp
    test/detail/output_string.cpp
    test/detail/engine_impl.cpp
    test/detail/intermediate_handler.cpp

    test/detail/typing/meta_check_context.cpp
    test/detail/typing/pos_map.cpp
    test/detail/typing/readable_field_traits.cpp
    test/detail/typing/row_traits.cpp

    test/impl/dt_to_string.cpp
    test/impl/ssl_context_with_default.cpp
    test/impl/variant_stream.cpp

    test/spotchecks/connection_use_after_move.cpp
    test/spotchecks/default_completion_tokens.cpp
    test/spotchecks/read_some_rows_static.cpp
    test/spotchecks/multifn.cpp
    test/spotchecks/execution_requests.cpp
    test/spotchecks/misc.cpp

    test/format_sql/formattable.cpp
    test/format_sql/basic_format_context.cpp
    test/format_sql/individual_value.cpp
    test/format_sql/formattable_ref.cpp
    test/format_sql/ranges.cpp
    test/format_sql/sequence.cpp
    test/format_sql/custom_formatter.cpp
    test/format_sql/format_strings.cpp
    test/format_sql/api.cpp

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
    test/is_fatal_error.cpp
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
    test/any_address.cpp
    test/any_connection.cpp
    test/pool_params.cpp
    test/connection_pool.cpp
    test/character_set.cpp
    test/escape_string.cpp
    test/constant_string_view.cpp
    test/pfr.cpp
    test/pipeline.cpp
    test/with_diagnostics.cpp
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

add_test(
    NAME boost_mysql_unittests
    COMMAND boost_mysql_unittests
)
