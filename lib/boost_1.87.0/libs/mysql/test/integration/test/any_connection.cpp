//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/engine_impl.hpp>

#include <boost/mysql/impl/internal/variant_stream.hpp>

#include <boost/asio/as_tuple.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/local/basic_endpoint.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/test/data/test_case.hpp>

#include <chrono>
#include <string>

#include "test_common/create_basic.hpp"
#include "test_common/create_diagnostics.hpp"
#include "test_common/network_result.hpp"
#include "test_common/poll_until.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/run_coro.hpp"
#include "test_integration/server_features.hpp"
#include "test_integration/spotchecks_helpers.hpp"

// Additional spotchecks for any_connection

using namespace boost::mysql;
using namespace boost::mysql::test;
namespace asio = boost::asio;
using boost::test_tools::per_element;

BOOST_AUTO_TEST_SUITE(test_any_connection)

// any_connection can be used with UNIX sockets
#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
BOOST_TEST_DECORATOR(*run_if(&server_features::unix_sockets))
BOOST_DATA_TEST_CASE(unix_sockets, network_functions_any::sync_and_async())
{
    // Setup
    netfn_fixture_any fix(sample);

    // Connect
    fix.connect(connect_params_builder().set_unix());
    BOOST_TEST(!fix.conn.uses_ssl());

    // We can prepare statements
    auto stmt = fix.net.prepare_statement(fix.conn, "SELECT ?, ?").get();
    BOOST_TEST(stmt.num_params() == 2u);

    // We can execute queries
    results r;
    fix.net.execute_query(fix.conn, "SELECT 'abc'", r).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, "abc"), per_element());

    // We can execute statements
    fix.net.execute_statement(fix.conn, stmt.bind(42, 100), r).validate_no_error();
    BOOST_TEST(r.rows() == makerows(2, 42, 100), per_element());

    // We can get errors
    fix.net.execute_query(fix.conn, "SELECT * FROM bad_table", r)
        .validate_error(
            common_server_errc::er_no_such_table,
            "Table 'boost_mysql_integtests.bad_table' doesn't exist"
        );

    // We can terminate the connection
    fix.net.close(fix.conn).validate_no_error();
}
#else
BOOST_DATA_TEST_CASE(unix_sockets_not_supported, network_functions_any::sync_and_async())
{
    // Setup
    netfn_fixture_any fix(sample);

    // Attempting to connect yields an error
    fix.net.connect(fix.conn, connect_params_builder().set_unix().build())
        .validate_error(asio::error::operation_not_supported);
}
#endif

// Backslash escapes
BOOST_FIXTURE_TEST_CASE(backslash_escapes, any_connection_fixture)
{
    // Backslash escapes enabled by default
    BOOST_TEST(conn.backslash_escapes());

    // Connect doesn't change the value
    conn.async_connect(connect_params_builder().disable_ssl().build(), as_netresult).validate_no_error();
    BOOST_TEST(conn.backslash_escapes());
    BOOST_TEST(conn.format_opts()->backslash_escapes);

    // Setting the SQL mode to NO_BACKSLASH_ESCAPES updates the value
    results r;
    conn.async_execute("SET sql_mode = 'NO_BACKSLASH_ESCAPES'", r, as_netresult).validate_no_error();
    BOOST_TEST(!conn.backslash_escapes());
    BOOST_TEST(!conn.format_opts()->backslash_escapes);

    // Executing a different statement doesn't change the value
    conn.async_execute("SELECT 1", r, as_netresult).validate_no_error();
    BOOST_TEST(!conn.backslash_escapes());
    BOOST_TEST(!conn.format_opts()->backslash_escapes);

    // Clearing the SQL mode updates the value
    conn.async_execute("SET sql_mode = ''", r, as_netresult).validate_no_error();
    BOOST_TEST(conn.backslash_escapes());
    BOOST_TEST(conn.format_opts()->backslash_escapes);

    // Reconnecting clears the value
    conn.async_execute("SET sql_mode = 'NO_BACKSLASH_ESCAPES'", r, as_netresult).validate_no_error();
    BOOST_TEST(!conn.backslash_escapes());
    BOOST_TEST(!conn.format_opts()->backslash_escapes);
    conn.async_connect(connect_params_builder().disable_ssl().build(), as_netresult).validate_no_error();
    BOOST_TEST(conn.backslash_escapes());
    BOOST_TEST(conn.format_opts()->backslash_escapes);
}

// Max buffer sizes
BOOST_AUTO_TEST_CASE(max_buffer_size)
{
    // Create the connection
    any_connection_params params;
    params.initial_buffer_size = 512u;
    params.max_buffer_size = 512u;
    any_connection_fixture fix(params);

    // Connect
    fix.conn.async_connect(connect_params_builder().disable_ssl().build(), as_netresult).validate_no_error();

    // Reading and writing almost 512 bytes works
    results r;
    auto q = format_sql(fix.conn.format_opts().value(), "SELECT {}", std::string(450, 'a'));
    fix.conn.async_execute(q, r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, std::string(450, 'a')), per_element());

    // Trying to write more than 512 bytes fails
    q = format_sql(fix.conn.format_opts().value(), "SELECT LENGTH({})", std::string(512, 'a'));
    fix.conn.async_execute(q, r, as_netresult).validate_error(client_errc::max_buffer_size_exceeded);

    // Trying to read more than 512 bytes fails
    fix.conn.async_execute("SELECT REPEAT('a', 512)", r, as_netresult)
        .validate_error(client_errc::max_buffer_size_exceeded);
}

BOOST_FIXTURE_TEST_CASE(default_max_buffer_size_success, any_connection_fixture)
{
    // Connect
    conn.async_connect(connect_params_builder().disable_ssl().build(), as_netresult).validate_no_error();

    // Reading almost max_buffer_size works
    execution_state st;
    conn.async_start_execution("SELECT 1, REPEAT('a', 0x3f00000)", st, as_netresult).validate_no_error();
    auto rws = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rws.at(0).at(1).as_string().size() == 0x3f00000u);
}

BOOST_FIXTURE_TEST_CASE(default_max_buffer_size_error, any_connection_fixture)
{
    // Connect
    conn.async_connect(connect_params_builder().disable_ssl().build(), as_netresult).validate_no_error();

    // Trying to read more than max_buffer_size bytes fails
    results r;
    conn.async_execute("SELECT 1, REPEAT('a', 0x4000000)", r, as_netresult)
        .validate_error(client_errc::max_buffer_size_exceeded);
}

BOOST_DATA_TEST_CASE(naggle_disabled, network_functions_any::sync_and_async())
{
    // Setup
    netfn_fixture_any fix(sample);

    // Connect
    fix.net.connect(fix.conn, connect_params_builder().disable_ssl().build()).validate_no_error();

    // Naggle's algorithm was disabled
    asio::ip::tcp::no_delay opt;
    static_cast<detail::engine_impl<detail::variant_stream>&>(detail::access::get_impl(fix.conn).get_engine())
        .stream()
        .socket()
        .get_option(opt);
    BOOST_TEST(opt.value() == true);
}

// Regression test: using a non-connected connection doesn't crash
BOOST_FIXTURE_TEST_CASE(using_non_connected_connection, any_connection_fixture)
{
    conn.async_ping(as_netresult).validate_any_error();
}

// Spotcheck: we can use cancel_after and other tokens
// that require initiations to have an associated executor
BOOST_FIXTURE_TEST_CASE(cancel_after, any_connection_fixture)
{
    // The token to use
    const auto token = asio::cancel_after(std::chrono::seconds(10), asio::deferred);

    // Connect
    conn.async_connect(connect_params_builder().build(), token)(as_netresult).validate_no_error();

    // Execute
    results result;
    conn.async_execute("SELECT 'abc'", result, token)(as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, "abc"), per_element());

    // Start execution
    execution_state st;
    conn.async_start_execution("SELECT 'abc'", st, token)(as_netresult).validate_no_error();
    auto rws = conn.async_read_some_rows(st, token)(as_netresult).get();
    BOOST_TEST(rws == makerows(1, "abc"), per_element());
    conn.async_read_resultset_head(st, token)(as_netresult).validate_no_error();

#ifdef BOOST_MYSQL_CXX14
    // Start execution (static, for read_some_rows)
    using tup_t = std::tuple<std::string>;
    static_execution_state<tup_t> st2;
    std::array<tup_t, 2> storage;
    conn.async_start_execution("SELECT 'abc'", st2, token)(as_netresult).validate_no_error();
    std::size_t sz = conn.async_read_some_rows(st2, boost::span<tup_t>(storage), token)(as_netresult).get();
    BOOST_TEST(sz == 1u);
#endif

    // Prepare & close statement
    auto stmt = conn.async_prepare_statement("SELECT ?", token)(as_netresult).get();
    conn.async_close_statement(stmt, token)(as_netresult).validate_no_error();

    // Reset connection & ping
    conn.async_reset_connection(token)(as_netresult).validate_no_error();
    conn.async_ping(token)(as_netresult).validate_no_error();

    // Close
    conn.async_close(token)(as_netresult).validate_no_error();
}

#ifdef BOOST_ASIO_HAS_CO_AWAIT

// Spotcheck: we can co_await async functions in any_connection,
// and this throws the right exception type
BOOST_FIXTURE_TEST_CASE(default_token, any_connection_fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        // Connect
        co_await conn.async_connect(connect_params_builder().build());

        // Success case
        results result;
        co_await conn.async_execute("SELECT 'abc'", result);
        BOOST_TEST(result.rows() == makerows(1, "abc"), per_element());

        // Error case
        BOOST_CHECK_EXCEPTION(
            co_await conn.async_execute("SELECT * FROM bad_table", result),
            error_with_diagnostics,
            [](const error_with_diagnostics& err) {
                BOOST_TEST(err.code() == common_server_errc::er_no_such_table);
                BOOST_TEST(
                    err.get_diagnostics() ==
                    create_server_diag("Table 'boost_mysql_integtests.bad_table' doesn't exist")
                );
                return true;
            }
        );

        // Returning a value works
        auto stmt = co_await conn.async_prepare_statement("SELECT ?");
        BOOST_TEST(stmt.valid());
    });
}

// The pattern co_await conn.fn(..., cancel_after(10s)) works
BOOST_FIXTURE_TEST_CASE(default_token_cancel_after, any_connection_fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        constexpr std::chrono::seconds timeout(10);

        // Connect
        co_await conn.async_connect(connect_params_builder().build(), asio::cancel_after(timeout));

        // Returning a value works
        auto stmt = co_await conn.async_prepare_statement("SELECT ?", asio::cancel_after(timeout));
        BOOST_TEST(stmt.valid());

        // Error case
        results result;
        BOOST_CHECK_EXCEPTION(
            co_await conn.async_execute("SELECT * FROM bad_table", result, asio::cancel_after(timeout)),
            error_with_diagnostics,
            [](const error_with_diagnostics& err) {
                BOOST_TEST(err.code() == common_server_errc::er_no_such_table);
                BOOST_TEST(
                    err.get_diagnostics() ==
                    create_server_diag("Table 'boost_mysql_integtests.bad_table' doesn't exist")
                );
                return true;
            }
        );
    });
}

// Using as_tuple as partial token works
BOOST_FIXTURE_TEST_CASE(default_token_as_tuple, any_connection_fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        // connect
        auto [ec] = co_await conn.async_connect(connect_params_builder().build(), asio::as_tuple);
        BOOST_TEST_REQUIRE(ec == error_code());

        // Returning a value works
        auto [ec2, stmt] = co_await conn.async_prepare_statement("SELECT ?", asio::as_tuple);
        BOOST_TEST_REQUIRE(ec2 == error_code());
        BOOST_TEST(stmt.valid());

        // Error case
        results result;
        auto [ec3] = co_await conn.async_execute("SELECT * FROM bad_table", result, asio::as_tuple);
        BOOST_TEST(ec3 == common_server_errc::er_no_such_table);
    });
}

// Using redirect_error as partial token works
BOOST_FIXTURE_TEST_CASE(default_token_redirect_error, any_connection_fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        // connect
        error_code ec;
        co_await conn.async_connect(connect_params_builder().build(), asio::redirect_error(ec));
        BOOST_TEST_REQUIRE(ec == error_code());

        // Returning a value works
        auto stmt = co_await conn.async_prepare_statement("SELECT ?", asio::redirect_error(ec));
        BOOST_TEST_REQUIRE(ec == error_code());
        BOOST_TEST(stmt.valid());

        // Error case
        results result;
        co_await conn.async_execute("SELECT * FROM bad_table", result, asio::redirect_error(ec));
        BOOST_TEST(ec == common_server_errc::er_no_such_table);
    });
}

// Spotcheck: immediate completions dispatched to the immediate executor
BOOST_FIXTURE_TEST_CASE(immediate_completions, any_connection_fixture)
{
    run_in_context(ctx, [this]() {
        // Setup
        connect();
        results r;

        // Prepare a statement
        auto stmt = conn.async_prepare_statement("SELECT 1", as_netresult).get();

        // Executing with the wrong number of params is an immediate error
        conn.async_execute(stmt.bind(0), r, as_netresult)
            .run()
            .validate_immediate(true)
            .validate_error(client_errc::wrong_num_params);
    });
}

#endif

BOOST_AUTO_TEST_SUITE_END()