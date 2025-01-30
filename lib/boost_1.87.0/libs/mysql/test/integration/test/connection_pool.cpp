//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/strand.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cstddef>
#include <exception>
#include <memory>
#include <stdexcept>
#include <utility>

#include "test_common/ci_server.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_common/tracker_executor.hpp"
#include "test_integration/run_coro.hpp"
#include "test_integration/server_features.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::test_tools::per_element;
namespace data = boost::unit_test::data;
namespace asio = boost::asio;

BOOST_AUTO_TEST_SUITE(test_connection_pool)

pool_params create_pool_params(std::size_t max_size = 151)
{
    pool_params res;
    res.server_address.emplace_host_and_port(get_hostname());
    res.username = integ_user;
    res.password = integ_passwd;
    res.database = integ_db;
    res.ssl = ssl_mode::disable;
    res.max_size = max_size;
    return res;
}

static void check_run(error_code ec)
{
    // Should complete successfully
    BOOST_TEST(ec == error_code());

    // Should never complete immediately
    BOOST_TEST(!is_initiation_function());
}

struct fixture : io_context_fixture
{
    // async_get_connection actually passes nullptr as diagnostics* to initiation
    // functions if no diagnostics is provided (unlike any_connection)
    diagnostics diag;
    results r;
};

// The pool and individual connections use the correct executors
BOOST_FIXTURE_TEST_CASE(connection_executor, fixture)
{
    // Create two different executors
    asio::any_io_executor pool_ex = asio::make_strand(ctx);
    asio::any_io_executor conn_ex = ctx.get_executor();
    BOOST_TEST((pool_ex != conn_ex));

    // Create and run the pool
    auto params = create_pool_params();
    params.connection_executor = conn_ex;
    connection_pool pool(pool_ex, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Check executors
    BOOST_TEST((pool.get_executor() == pool_ex));
    BOOST_TEST((conn->get_executor() == conn_ex));

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

BOOST_FIXTURE_TEST_CASE(pool_executors_thread_safe, fixture)
{
    // Create and run the pool
    auto params = create_pool_params();
    params.thread_safe = true;
    connection_pool pool(ctx, create_pool_params());
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Check executors. The internal strand is never exposed,
    // and doesn't get propagated to connections
    BOOST_TEST((pool.get_executor() == ctx.get_executor()));
    BOOST_TEST((conn->get_executor() == ctx.get_executor()));

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

BOOST_DATA_TEST_CASE_F(fixture, return_connection_with_reset, data::make({false, true}))
{
    // Create a pool with max_size 1, so the same connection gets always returned
    auto params = create_pool_params(1);
    params.thread_safe = sample;
    connection_pool pool(ctx, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Alter session state
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_execute("SET @myvar = 'abc'", r, as_netresult).validate_no_error();

    // Return the connection
    conn = pooled_connection();

    // Get the same connection again
    conn = pool.async_get_connection(diag, as_netresult).get();

    // The same connection is returned, but session state has been cleared
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_execute("SELECT @myvar", r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, nullptr), per_element());

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

BOOST_DATA_TEST_CASE_F(fixture, return_connection_without_reset, data::make({false, true}))
{
    // Create a connection pool with max_size 1, so the same connection gets always returned
    auto params = create_pool_params(1);
    params.thread_safe = sample;
    connection_pool pool(ctx, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Alter session state
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_execute("SET @myvar = 'abc'", r, as_netresult).validate_no_error();

    // Return the connection
    conn.return_without_reset();
    BOOST_TEST(!conn.valid());

    // Get the same connection again
    conn = pool.async_get_connection(diag, as_netresult).get();

    // The same connection is returned, and no reset has been issued
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_execute("SELECT @myvar", r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, "abc"), per_element());

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

// pooled_connection destructor is equivalent to return_connection with reset
BOOST_FIXTURE_TEST_CASE(pooled_connection_destructor, fixture)
{
    // Create a connection pool with max_size 1, so the same connection gets always returned
    connection_pool pool(ctx, create_pool_params(1));
    auto run_result = pool.async_run(as_netresult);

    {
        // Get a connection
        auto conn = pool.async_get_connection(diag, as_netresult).get();

        // Alter session state
        BOOST_TEST_REQUIRE(conn.valid());
        conn->async_execute("SET @myvar = 'abc'", r, as_netresult).validate_no_error();
    }

    // Get the same connection again
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // The same connection is returned, but session state has been cleared
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_execute("SELECT @myvar", r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, nullptr), per_element());

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

// Pooled connections use utf8mb4
static void validate_charset(any_connection& conn)
{
    // The connection knows its using utf8mb4
    BOOST_TEST(conn.current_character_set()->name == "utf8mb4");
    BOOST_TEST(conn.format_opts()->charset.name == "utf8mb4");

    // The connection is actually using utf8mb4
    results r;
    conn.async_execute(
            "SELECT @@character_set_client, @@character_set_connection, @@character_set_results",
            r,
            as_netresult
    )
        .validate_no_error();
    BOOST_TEST(r.rows() == makerows(3, "utf8mb4", "utf8mb4", "utf8mb4"), per_element());
}

BOOST_FIXTURE_TEST_CASE(charset, fixture)
{
    // Create and run the pool
    connection_pool pool(ctx, create_pool_params(1));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();
    validate_charset(conn.get());

    // Return the connection and retrieve it again
    conn = pooled_connection();
    conn = pool.async_get_connection(diag, as_netresult).get();
    validate_charset(conn.get());

    // Return the connection without reset and retrieve it again
    conn.return_without_reset();
    conn = pool.async_get_connection(diag, as_netresult).get();
    validate_charset(conn.get());

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

BOOST_FIXTURE_TEST_CASE(connections_created_if_required, fixture)
{
    connection_pool pool(ctx, create_pool_params());
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn1 = pool.async_get_connection(diag, as_netresult).get();

    // Check that it works
    BOOST_TEST_REQUIRE(conn1.valid());
    conn1->async_execute("SET @myvar = '1'", r, as_netresult).validate_no_error();

    // Get another connection. This will create a new one, since the first one is in use
    auto conn2 = pool.async_get_connection(diag, as_netresult).get();

    // Check that it works
    BOOST_TEST_REQUIRE(conn1.valid());
    conn2->async_execute("SET @myvar = '2'", r, as_netresult).validate_no_error();

    // They are different connections
    conn1->async_execute("SELECT @myvar", r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, "1"), per_element());
    conn2->async_execute("SELECT @myvar", r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, "2"), per_element());

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

BOOST_FIXTURE_TEST_CASE(connection_upper_limit, fixture)
{
    connection_pool pool(ctx, create_pool_params(1));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Getting another connection will block until one is returned.
    // Since we won't return the one we have, the function time outs
    pool.async_get_connection(diag, asio::cancel_after(std::chrono::milliseconds(1), asio::deferred))(
            as_netresult
    )
        .validate_error(client_errc::no_connection_available);

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

// If a connection is requested before calling run, we wait
BOOST_DATA_TEST_CASE_F(fixture, get_connection_before_run, data::make({false, true}))
{
    auto params = create_pool_params(1);
    params.thread_safe = sample;
    connection_pool pool(ctx, std::move(params));

    // Get a connection before calling run
    auto getconn_result = pool.async_get_connection(diag, as_netresult);

    // Call run
    auto run_result = pool.async_run(as_netresult);

    // Success
    auto conn = std::move(getconn_result).get();
    conn->async_ping(as_netresult).validate_no_error();

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

BOOST_FIXTURE_TEST_CASE(cancel_run, fixture)
{
    // Construct a pool and run it
    connection_pool pool(ctx, create_pool_params());
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Cancel. This will make run() return
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();

    // Cancel again does nothing
    pool.cancel();
}

// If the pool is cancelled before calling run, cancel still has effect
BOOST_FIXTURE_TEST_CASE(cancel_before_run, fixture)
{
    // Create a pool
    connection_pool pool(ctx, create_pool_params());

    // Cancel
    pool.cancel();

    // Run returns immediately
    pool.async_run(as_netresult).validate_no_error_nodiag();
}

BOOST_DATA_TEST_CASE_F(fixture, cancel_get_connection, data::make({false, true}))
{
    // Construct a pool and run it
    auto params = create_pool_params(1);
    params.thread_safe = sample;
    connection_pool pool(ctx, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Try to get a new one. This will not complete, since there is no room for more connections
    diagnostics diag2;
    auto getconn_result = pool.async_get_connection(diag2, as_netresult);

    // Cancel. This will make run and get_connection return
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
    std::move(getconn_result).validate_error(client_errc::pool_cancelled);

    // Calling get_connection after cancel will error
    pool.async_get_connection(diag, as_netresult).validate_error(client_errc::pool_cancelled);
}

// Connection pool's destructor cancels the pool
BOOST_FIXTURE_TEST_CASE(destructor_cancel, fixture)
{
    // Construct a pool and run it
    std::unique_ptr<connection_pool> pool{new connection_pool(ctx, create_pool_params(1))};
    auto run_result = pool->async_run(as_netresult);

    // Try to get 2 connections. The 2nd one blocks
    auto conn = pool->async_get_connection(diag, as_netresult).get();
    auto getconn_result = pool->async_get_connection(diag, as_netresult);

    // Destroy the pool
    pool.reset();

    // Run returns and the connection request is cancelled
    std::move(run_result).validate_no_error_nodiag();
    std::move(getconn_result).validate_error(client_errc::pool_cancelled);
}

// Having a valid pooled_connection alive extends the pool's lifetime
BOOST_DATA_TEST_CASE_F(fixture, pooled_connection_extends_pool_lifetime, data::make({false, true}))
{
    auto params = create_pool_params();
    params.thread_safe = sample;
    std::unique_ptr<connection_pool> pool(new connection_pool(ctx, std::move(params)));

    // Run the pool
    auto run_result = pool->async_run(as_netresult);

    // Get a connection
    auto conn = pool->async_get_connection(diag, as_netresult).get();

    // Cancel and destroy
    pool->cancel();
    pool.reset();

    // Wait for run to exit, since run extends lifetime, too
    std::move(run_result).validate_no_error_nodiag();

    // The connection we got can still be used and returned
    // In thread-safe mode, strand dispatching doesn't cause lifetime problems
    conn->async_ping(as_netresult).validate_no_error();
    conn.return_without_reset();
}

// Having a packaged async_get_connection op extends lifetime
BOOST_FIXTURE_TEST_CASE(async_get_connection_initation_extends_pool_lifetime, fixture)
{
    std::unique_ptr<connection_pool> pool(new connection_pool(ctx, create_pool_params()));

    // Create a packaged op
    auto op = pool->async_get_connection(diag, boost::asio::deferred);

    // Destroy the pool
    pool.reset();

    // We can run the operation without crashing, since it extends lifetime
    std::move(op)(asio::cancel_after(std::chrono::nanoseconds(1), as_netresult))
        .validate_error(client_errc::pool_not_running);
}

// In thread-safe mode, cancel() is dispatched to the strand, and doesn't cause lifetime issues
BOOST_FIXTURE_TEST_CASE(cancel_extends_pool_lifetime, fixture)
{
    auto params = create_pool_params();
    params.thread_safe = true;
    std::unique_ptr<connection_pool> pool(new connection_pool(ctx, std::move(params)));

    // Cancel
    pool->cancel();

    // Destroy the pool
    pool.reset();

    // Dispatch any pending handler. We didn't crash
    ctx.poll();
}

// Spotcheck: the overload without diagnostics work
BOOST_FIXTURE_TEST_CASE(get_connection_no_diag, fixture)
{
    connection_pool pool(ctx, create_pool_params());
    auto run_result = pool.async_run(as_netresult);

    auto conn = pool.async_get_connection(as_netresult).get_nodiag();
    conn->async_ping(as_netresult).validate_no_error();

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
// Spotcheck: pool works with unix sockets, too
BOOST_TEST_DECORATOR(*run_if(&server_features::unix_sockets))
BOOST_FIXTURE_TEST_CASE(unix_sockets, fixture)
{
    // Create and run the pool
    auto params = create_pool_params();
    params.server_address.emplace_unix_path(default_unix_path);
    connection_pool pool(ctx, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Verify that works
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_ping(as_netresult).validate_no_error();

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}
#endif

// Spotcheck: pool works with TLS
BOOST_FIXTURE_TEST_CASE(ssl, fixture)
{
    // Create and run the pool
    auto params = create_pool_params();
    params.ssl = ssl_mode::require;
    connection_pool pool(ctx, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Verify that works
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_ping(as_netresult).validate_no_error();

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

// Spotcheck: custom ctor params (SSL context and buffer size) can be passed to the connection pool
BOOST_FIXTURE_TEST_CASE(custom_ctor_params, fixture)
{
    // Create and run the pool
    auto params = create_pool_params();
    params.ssl = ssl_mode::require;
    params.ssl_ctx.emplace(asio::ssl::context::sslv23_client);
    params.initial_buffer_size = 16u;
    connection_pool pool(ctx, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Verify that works
    BOOST_TEST_REQUIRE(conn.valid());
    conn->async_ping(as_netresult).validate_no_error();

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

// Spotcheck: the pool can work with zero timeouts
BOOST_FIXTURE_TEST_CASE(zero_timeuts, fixture)
{
    // Create and run the pool
    auto params = create_pool_params();
    params.max_size = 1u;  // so we force a reset
    params.connect_timeout = std::chrono::seconds(0);
    params.ping_timeout = std::chrono::seconds(0);
    params.ping_interval = std::chrono::seconds(0);
    connection_pool pool(ctx, std::move(params));
    auto run_result = pool.async_run(as_netresult);

    // Get a connection
    auto conn = pool.async_get_connection(diag, as_netresult).get();
    conn->async_ping(as_netresult).validate_no_error();

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

// Spotcheck: we can use completion tokens that require
// initiations to have a bound executor, like cancel_after
// This also tests that running ops with a connected cancel slot
// without triggering cancellation doesn't crash
BOOST_FIXTURE_TEST_CASE(cancel_after, fixture)
{
    constexpr std::chrono::seconds timeout(10);

    connection_pool pool(ctx, create_pool_params());
    pool.async_run(asio::cancel_after(timeout, check_run));

    // Get a connection
    auto conn = pool.async_get_connection(diag, asio::cancel_after(timeout, asio::deferred))(as_netresult)
                    .get();
    conn->async_ping(as_netresult).validate_no_error();

    // Cleanup the pool
    pool.cancel();
}

// Spotcheck: per-operation cancellation works with async_run
BOOST_FIXTURE_TEST_CASE(async_run_per_operation_cancellation, fixture)
{
    connection_pool pool(ctx, create_pool_params());
    pool.async_run(asio::cancel_after(std::chrono::microseconds(1), asio::deferred))(as_netresult)
        .validate_no_error_nodiag();
    pool.async_get_connection(diag, as_netresult).validate_error(client_errc::pool_cancelled);
}

// Spotcheck: per-operation cancellation works with async_get_connection
BOOST_FIXTURE_TEST_CASE(async_get_connection_per_operation_cancellation, fixture)
{
    // Create and run the pool
    connection_pool pool(ctx, create_pool_params(1));
    auto run_result = pool.async_run(as_netresult);

    // Get the only connection the pool has
    auto conn = pool.async_get_connection(diag, as_netresult).get();

    // Getting another connection times out
    pool.async_get_connection(diag, asio::cancel_after(std::chrono::microseconds(1), asio::deferred))(
            as_netresult
    )
        .validate_error(client_errc::no_connection_available);

    // Cleanup the pool
    pool.cancel();
    std::move(run_result).validate_no_error_nodiag();
}

#ifdef BOOST_ASIO_HAS_CO_AWAIT

// Spotcheck: we can co_await async functions in any_connection,
// and this throws the right exception type
BOOST_FIXTURE_TEST_CASE(default_token, fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        connection_pool pool(ctx, create_pool_params());

        // Run can be used without a token. Defaults to with_diagnostics(deferred)
        auto run_op = pool.async_run();

        // Error case (pool not running)
        BOOST_CHECK_EXCEPTION(
            co_await pool.async_get_connection(asio::cancel_after(std::chrono::nanoseconds(1))),
            error_with_diagnostics,
            [](const error_with_diagnostics& err) {
                BOOST_TEST(err.code() == client_errc::pool_not_running);
                BOOST_TEST(err.get_diagnostics() == diagnostics());
                return true;
            }
        );

        // Run the pool
        std::move(run_op)([](std::exception_ptr exc) {
            if (exc)
                std::rethrow_exception(exc);
        });

        // Success case
        auto conn = co_await pool.async_get_connection();
        co_await conn->async_ping();

        // Finish
        pool.cancel();
    });
}

// cancel_after can be used as a partial token with async_run and async_get_connection.
BOOST_FIXTURE_TEST_CASE(cancel_after_partial_token, fixture)
{
    run_coro(ctx, [&]() -> asio::awaitable<void> {
        connection_pool pool(ctx, create_pool_params(1));

        // Run can be used with cancel_after
        asio::co_spawn(
            ctx,
            [&]() -> asio::awaitable<void> {
                co_await pool.async_run(asio::cancel_after(std::chrono::seconds(1)));
            },
            [](std::exception_ptr exc) {
                if (exc)
                    std::rethrow_exception(exc);
            }
        );

        // Success case
        auto conn = co_await pool.async_get_connection(asio::cancel_after(std::chrono::seconds(1)));
        co_await conn->async_ping();

        // Error case (operation cancelled)
        BOOST_CHECK_EXCEPTION(
            co_await pool.async_get_connection(asio::cancel_after(std::chrono::nanoseconds(1))),
            error_with_diagnostics,
            [](const error_with_diagnostics& err) {
                BOOST_TEST(err.code() == client_errc::no_connection_available);
                BOOST_TEST(err.get_diagnostics() == diagnostics());
                return true;
            }
        );

        // Finish
        pool.cancel();
    });
}

#endif

// Spotcheck: constructing a connection_pool with invalid params throws
BOOST_AUTO_TEST_CASE(invalid_params)
{
    asio::io_context ctx;
    pool_params params;
    params.connect_timeout = std::chrono::seconds(-1);

    BOOST_CHECK_EXCEPTION(
        connection_pool(ctx, std::move(params)),
        std::invalid_argument,
        [](const std::invalid_argument& exc) {
            return exc.what() == string_view("pool_params::connect_timeout must not be negative");
        }
    );
}

BOOST_AUTO_TEST_SUITE_END()
