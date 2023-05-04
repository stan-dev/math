//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/protocol/resultset_encoding.hpp>

#include <boost/asio/use_awaitable.hpp>
#include <boost/test/unit_test.hpp>

#include <tuple>

#include "assert_buffer_equals.hpp"
#include "create_execution_state.hpp"
#include "create_message.hpp"
#include "create_statement.hpp"
#include "printing.hpp"
#include "run_coroutine.hpp"
#include "test_common.hpp"
#include "test_connection.hpp"
#include "unit_netfun_maker.hpp"

using boost::mysql::blob;
using boost::mysql::client_errc;
using boost::mysql::column_type;
using boost::mysql::error_code;
using boost::mysql::results;
using boost::mysql::statement;
using boost::mysql::detail::connection_access;
using boost::mysql::detail::protocol_field_type;
using boost::mysql::detail::results_access;
using boost::mysql::detail::resultset_encoding;
using namespace boost::mysql::test;

namespace {

// Machinery to treat iterator and tuple overloads the same
using netfun_maker = netfun_maker_mem<
    void,
    test_connection,
    const statement&,
    const std::tuple<const char*, std::nullptr_t>&,
    results&>;

struct
{
    netfun_maker::signature execute_statement;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&test_connection::execute_statement),             "sync_errc"      },
    {netfun_maker::sync_exc(&test_connection::execute_statement),              "sync_exc"       },
    {netfun_maker::async_errinfo(&test_connection::async_execute_statement),   "async_errinfo"  },
    {netfun_maker::async_noerrinfo(&test_connection::async_execute_statement), "async_noerrinfo"},
};

// Verify that we reset the results object
results create_initial_results()
{
    results res;
    results_access::get_rows(res) = makerows(1, 42, "abc");
    results_access::get_state(res
    ) = create_execution_state(resultset_encoding::text, {protocol_field_type::geometry}, 4);
    return res;
}

BOOST_AUTO_TEST_SUITE(test_execute_statement)

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto result = create_initial_results();
            test_connection conn;
            auto stmt = create_statement(2);
            conn.stream().add_message(create_ok_packet_message(1, 2, 3, 4, 5, "info"));

            // Call the function
            fns.execute_statement(conn, stmt, std::make_tuple("test", nullptr), result).validate_no_error();

            // Verify the message we sent
            std::uint8_t expected_message[] = {
                0x015, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
                0x00,  0x02, 0x01, 0xfe, 0x00, 0x06, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74,
            };
            BOOST_MYSQL_ASSERT_BLOB_EQUALS(conn.stream().bytes_written(), expected_message);

            // Verify the results
            BOOST_TEST(result.meta().size() == 0u);
            BOOST_TEST(result.affected_rows() == 2u);
            BOOST_TEST(result.last_insert_id() == 3u);
            BOOST_TEST(result.warning_count() == 5u);
            BOOST_TEST(result.info() == "info");
        }
    }
}

BOOST_AUTO_TEST_CASE(error_start_statement_execution)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto result = create_initial_results();
            test_connection conn;
            auto stmt = create_statement(2);
            conn.stream().set_fail_count(fail_count(0, client_errc::server_unsupported));

            // Call the function
            fns.execute_statement(conn, stmt, std::make_tuple("abc", nullptr), result)
                .validate_error_exact(client_errc::server_unsupported);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_read_all_rows)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            auto result = create_initial_results();
            test_connection conn;
            auto stmt = create_statement(2);
            conn.stream().add_message(create_message(1, {0x01}));  // Response OK, 1 metadata packet
            conn.stream().add_message(create_coldef_message(2, protocol_field_type::geometry));
            conn.stream().set_fail_count(fail_count(4, client_errc::server_unsupported));

            // Call the function
            fns.execute_statement(conn, stmt, std::make_tuple("abc", nullptr), result)
                .validate_error_exact(client_errc::server_unsupported);

            // Ensure we successfully ran the start_query
            BOOST_TEST_REQUIRE(results_access::get_state(result).meta().size() == 1u);
            BOOST_TEST(results_access::get_state(result).meta()[0].type() == column_type::geometry);
        }
    }
}

// Verify that we correctly perform a decay-copy of the parameters and the
// statement handle, relevant for deferred tokens
#ifdef BOOST_ASIO_HAS_CO_AWAIT
BOOST_AUTO_TEST_SUITE(deferred_lifetimes)
struct fixture
{
    results result{create_initial_results()};
    test_connection conn;
    statement stmt{create_statement(2)};

    static constexpr std::uint8_t expected_msg[]{
        0x1d, 0x00, 0x00, 0x00, 0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0xfe,
        0x00, 0x08, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74, 0x2a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    };

    fixture() { conn.stream().add_message(create_ok_packet_message(1)); }
};

BOOST_AUTO_TEST_CASE(params_rvalue)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        auto aw = fix.conn.async_execute_statement(
            fix.stmt,
            std::make_tuple(std::string("test"), 42),
            fix.result,
            boost::asio::use_awaitable
        );
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), fix.expected_msg);
        BOOST_TEST(fix.result.rows().size() == 0u);
    });
}

BOOST_AUTO_TEST_CASE(params_lvalue)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        auto tup = std::make_tuple(std::string("test"), 42);
        auto aw = fix.conn.async_execute_statement(fix.stmt, tup, fix.result, boost::asio::use_awaitable);
        tup = std::make_tuple(std::string("other"), 90);
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), fix.expected_msg);
        BOOST_TEST(fix.result.rows().size() == 0u);
    });
}

BOOST_AUTO_TEST_CASE(params_const_lvalue)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        const auto tup = std::make_tuple(std::string("test"), 42);
        auto aw = fix.conn.async_execute_statement(fix.stmt, tup, fix.result, boost::asio::use_awaitable);
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), fix.expected_msg);
        BOOST_TEST(fix.result.rows().size() == 0u);
    });
}

// Verify that we don't require the passed-in statement to be alive
BOOST_AUTO_TEST_CASE(statement_handle)
{
    run_coroutine([]() -> boost::asio::awaitable<void> {
        fixture fix;

        // Deferred op
        auto aw = fix.conn.async_execute_statement(
            statement(fix.stmt),
            std::make_tuple(std::string("test"), 42),
            fix.result,
            boost::asio::use_awaitable
        );
        co_await std::move(aw);

        // verify that the op had the intended effects
        BOOST_MYSQL_ASSERT_BLOB_EQUALS(fix.conn.stream().bytes_written(), fix.expected_msg);
        BOOST_TEST(fix.result.rows().size() == 0u);
    });
}

BOOST_AUTO_TEST_SUITE_END()
#endif

BOOST_AUTO_TEST_SUITE_END()

}  // namespace