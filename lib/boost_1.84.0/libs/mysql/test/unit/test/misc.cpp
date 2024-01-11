//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/mysql_collations.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/asio/use_awaitable.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <memory>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/create_statement.hpp"
#include "test_unit/run_coroutine.hpp"
#include "test_unit/test_stream.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_misc)

using test_connection = connection<test_stream>;

// Make sure async_execute() and friends don't cause side
// effects in the initiation
#ifdef BOOST_ASIO_HAS_CO_AWAIT
BOOST_AUTO_TEST_CASE(async_execute_side_effects_in_initiation)
{
    test_connection conn;
    results result1, result2;

    // Resultsets will be complete as soon as a message is read
    conn.stream()
        .add_bytes(create_ok_frame(1, ok_builder().affected_rows(2).build()))
        .add_bytes(create_ok_frame(1, ok_builder().affected_rows(1).build()));

    // Launch coroutine and wait for completion
    run_coroutine(conn.get_executor(), [&]() -> boost::asio::awaitable<void> {
        // Call both queries but don't wait on them yet, so they don't initiate
        auto aw1 = conn.async_query("Q1", result1, boost::asio::use_awaitable);
        auto aw2 = conn.async_query("Q2", result2, boost::asio::use_awaitable);

        // Run them in reverse order
        co_await std::move(aw2);
        co_await std::move(aw1);
    });

    // Check that we wrote Q2's message first, then Q1's
    auto expected = concat_copy(
        create_frame(0, {0x03, 0x51, 0x32}),  // query request Q2
        create_frame(0, {0x03, 0x51, 0x31})   // query request Q1
    );
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), expected);

    // Check that the results got the right ok_packets
    BOOST_TEST(result2.affected_rows() == 2u);
    BOOST_TEST(result1.affected_rows() == 1u);
}
#endif  // BOOST_ASIO_HAS_CO_AWAIT

// spotcheck for the dynamic interface
// Verifies that execute (dynamic interface) works when rows come in separate batches
// This is testing the interaction between the network algorithm and results
BOOST_AUTO_TEST_CASE(execute_multiple_batches)
{
    // Setup
    test_connection conn;
    results result;

    // Message sequence (each on its own read)
    conn.stream()
        .add_bytes(create_frame(1, {0x02}))  // OK, 2 columns
        .add_break()
        .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))  // meta
        .add_break()
        .add_bytes(create_coldef_frame(
            3,
            meta_builder().type(column_type::blob).collation_id(mysql_collations::binary).build_coldef()
        ))  // meta
        .add_break()
        .add_bytes(create_text_row_message(4, "abcd", makebv("\0\1\0")))  // row 1
        .add_break()
        .add_bytes(create_text_row_message(5, "defghi", makebv("\3\4\3\0")))  // row 2
        .add_break()
        .add_bytes(create_eof_frame(6, ok_builder().affected_rows(10u).info("1st").more_results(true).build())
        )
        .add_break()
        .add_bytes(create_ok_frame(7, ok_builder().affected_rows(20u).info("2nd").more_results(true).build()))
        .add_break()
        .add_bytes(create_frame(8, {0x01}))  // OK, 1 metadata
        .add_break()
        .add_bytes(create_coldef_frame(9, meta_builder().type(column_type::varchar).build_coldef()))  // meta
        .add_break()
        .add_bytes(create_text_row_message(10, "ab"))  // row 1
        .add_break()
        .add_bytes(create_eof_frame(11, ok_builder().affected_rows(30u).info("3rd").build()));

    // Call the function
    conn.execute("abc", result);

    // We've written the query request
    auto expected_msg = create_frame(0, {0x03, 0x61, 0x62, 0x63});  // ASCII "abc" (plus length)
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), expected_msg);

    // We've populated the results
    BOOST_TEST_REQUIRE(result.size() == 3u);
    BOOST_TEST(result[0].affected_rows() == 10u);
    BOOST_TEST(result[0].info() == "1st");
    BOOST_TEST(result[0].rows() == makerows(2, "abcd", makebv("\0\1\0"), "defghi", makebv("\3\4\3\0")));
    BOOST_TEST(result[1].affected_rows() == 20u);
    BOOST_TEST(result[1].info() == "2nd");
    BOOST_TEST(result[1].rows() == rows());
    BOOST_TEST(result[2].affected_rows() == 30u);
    BOOST_TEST(result[2].info() == "3rd");
    BOOST_TEST(result[2].rows() == makerows(1, "ab"));
}

// Regression check: execute statement with iterator range with a reference type that is convertible to
// field_view, but not equal to field_view
BOOST_AUTO_TEST_CASE(execute_stmt_iterator_reference_not_field_view)
{
    results result;
    auto stmt = statement_builder().id(1).num_params(2).build();
    test_connection conn;
    conn.stream().add_bytes(create_ok_frame(1, ok_builder().affected_rows(50).info("1st").build()));

    // Call the function
    std::vector<field> fields{field_view("test"), field_view()};
    conn.execute(stmt.bind(fields.begin(), fields.end()), result);

    // Verify the message we sent
    constexpr std::uint8_t expected_msg[] = {
        0x15, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
        0x00, 0x02, 0x01, 0xfe, 0x00, 0x06, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74,
    };
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), expected_msg);

    // Verify the results
    BOOST_TEST_REQUIRE(result.size() == 1u);
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 50u);
    BOOST_TEST(result.info() == "1st");
}

#ifdef BOOST_ASIO_HAS_CO_AWAIT

// The serialized form of executing a statement with ID=1, params=("test", nullptr)
constexpr std::uint8_t execute_stmt_msg[] = {
    0x15, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
    0x00, 0x02, 0x01, 0xfe, 0x00, 0x06, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74,
};

// Verify that we correctly perform a decay-copy of the execution request
// in async_execute(), relevant for deferred tokens
BOOST_AUTO_TEST_CASE(async_execute_deferred_lifetimes_rvalues)
{
    test_connection conn;

    run_coroutine(conn.get_executor(), [&]() -> boost::asio::awaitable<void> {
        results result;
        conn.stream().add_bytes(create_ok_frame(1, ok_builder().info("1st").build()));

        // Deferred op. Execution request is a temporary
        auto aw = conn.async_execute(
            statement_builder().id(1).num_params(2).build().bind(std::string("test"), nullptr),
            result,
            boost::asio::use_awaitable
        );
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), execute_stmt_msg);
        BOOST_TEST(result.info() == "1st");
    });
}

BOOST_AUTO_TEST_CASE(async_execute_deferred_lifetimes_lvalues)
{
    test_connection conn;

    run_coroutine(conn.get_executor(), [&]() -> boost::asio::awaitable<void> {
        results result;

        // Create a bound statement on the heap. This helps tooling detect memory errors
        using bound_stmt_t = bound_statement_tuple<std::tuple<std::string, std::nullptr_t>>;
        auto stmt = statement_builder().id(1).num_params(2).build();
        std::unique_ptr<bound_stmt_t> stmt_ptr{new bound_stmt_t{stmt.bind(std::string("test"), nullptr)}};

        // Messages
        conn.stream().add_bytes(create_ok_frame(1, ok_builder().info("1st").build()));

        // Deferred op
        auto aw = conn.async_execute(*stmt_ptr, result, boost::asio::use_awaitable);

        // Free the statement
        stmt_ptr.reset();

        // Actually run the op
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), execute_stmt_msg);
        BOOST_TEST(result.info() == "1st");
    });
}

// Verify that we correctly perform a decay-copy of the parameters and the
// statement handle for async_start_execution(), relevant for deferred tokens
BOOST_AUTO_TEST_CASE(async_start_execution_deferred_lifetimes_rvalues)
{
    test_connection conn;

    run_coroutine(conn.get_executor(), [&]() -> boost::asio::awaitable<void> {
        execution_state st;
        conn.stream().add_bytes(create_ok_frame(1, ok_builder().info("1st").build()));

        // Deferred op. Execution request is a temporary
        auto aw = conn.async_start_execution(
            statement_builder().id(1).num_params(2).build().bind(std::string("test"), nullptr),
            st,
            boost::asio::use_awaitable
        );
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), execute_stmt_msg);
        BOOST_TEST(st.info() == "1st");
    });
}

BOOST_AUTO_TEST_CASE(deferred_lifetimes_lvalues)
{
    test_connection conn;

    run_coroutine(conn.get_executor(), [&]() -> boost::asio::awaitable<void> {
        execution_state st;
        conn.stream().add_bytes(create_ok_frame(1, ok_builder().info("1st").build()));

        // Create a bound statement on the heap. This helps tooling detect memory errors
        using bound_stmt_t = bound_statement_tuple<std::tuple<std::string, std::nullptr_t>>;
        auto stmt = statement_builder().id(1).num_params(2).build();
        std::unique_ptr<bound_stmt_t> stmt_ptr{new bound_stmt_t{stmt.bind(std::string("test"), nullptr)}};

        // Deferred op
        auto aw = conn.async_start_execution(*stmt_ptr, st, boost::asio::use_awaitable);

        // Free the statement
        stmt_ptr.reset();

        // Actually run the op
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), execute_stmt_msg);
        BOOST_TEST(st.info() == "1st");
    });
}

// Verify that async_close_statement doesn't require the passed-in statement to be alive. Only
// relevant for deferred tokens.
BOOST_AUTO_TEST_CASE(async_close_statement_handle_deferred_tokens)
{
    test_connection conn;

    run_coroutine(conn.get_executor(), [&]() -> boost::asio::awaitable<void> {
        auto stmt = statement_builder().id(3).build();

        // Deferred op
        auto aw = conn.async_close_statement(stmt, boost::asio::use_awaitable);

        // Invalidate the original variable
        stmt = statement_builder().id(42).build();

        // Run the operation
        co_await std::move(aw);

        // verify that the op had the intended effects
        const auto expected_message = create_frame(0, {0x19, 0x03, 0x00, 0x00, 0x00});
        BOOST_MYSQL_ASSERT_BUFFER_EQUALS(conn.stream().bytes_written(), expected_message);
    });
}
#endif

BOOST_AUTO_TEST_SUITE_END()
