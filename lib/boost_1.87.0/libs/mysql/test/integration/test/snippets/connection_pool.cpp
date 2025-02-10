//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/with_diagnostics.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/use_future.hpp>
#include <boost/test/unit_test.hpp>

#include <chrono>

#include "test_common/ci_server.hpp"
#include "test_integration/run_coro.hpp"
#include "test_integration/snippets/credentials.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

namespace {

#ifdef BOOST_ASIO_HAS_CO_AWAIT
//[connection_pool_get_connection
// Use connection pools for functions that will be called
// repeatedly during the application lifetime.
// An HTTP server handler function is a good candidate.
boost::asio::awaitable<std::int64_t> get_num_employees(boost::mysql::connection_pool& pool)
{
    // Get a fresh connection from the pool.
    // pooled_connection is a proxy to an any_connection object.
    boost::mysql::pooled_connection conn = co_await pool.async_get_connection();

    // Use pooled_connection::operator-> to access the underlying any_connection.
    // Let's use the connection
    results result;
    co_await conn->async_execute("SELECT COUNT(*) FROM employee", result);
    co_return result.rows().at(0).at(0).as_int64();

    // When conn is destroyed, the connection is returned to the pool
}
//]

boost::asio::awaitable<void> return_without_reset(boost::mysql::connection_pool& pool)
{
    //[connection_pool_return_without_reset
    // Get a connection from the pool
    boost::mysql::pooled_connection conn = co_await pool.async_get_connection();

    // Use the connection in a way that doesn't mutate session state.
    // We're not setting variables, preparing statements or starting transactions,
    // so it's safe to skip reset
    boost::mysql::results result;
    co_await conn->async_execute("SELECT COUNT(*) FROM employee", result);

    // Explicitly return the connection to the pool, skipping reset
    conn.return_without_reset();
    //]
}

boost::asio::awaitable<void> apply_timeout(boost::mysql::connection_pool& pool)
{
    //[connection_pool_apply_timeout
    // Get a connection from the pool, but don't wait more than 5 seconds
    auto conn = co_await pool.async_get_connection(boost::asio::cancel_after(std::chrono::seconds(5)));
    //]

    conn.return_without_reset();
}
#endif

//[connection_pool_sync
// Wraps a connection_pool and offers a sync interface.
// sync_pool is thread-safe
class sync_pool
{
    // A thread pool with a single thread. This is used to
    // run the connection pool. The thread is automatically
    // joined when sync_pool is destroyed.
    boost::asio::thread_pool thread_pool_{1};

    // The async connection pool
    boost::mysql::connection_pool conn_pool_;

public:
    // Constructor: constructs the connection_pool object from
    // the single-thread pool and calls async_run.
    // The pool has a single thread, which creates an implicit strand.
    // There is no need to use pool_params::thread_safe
    sync_pool(boost::mysql::pool_params params) : conn_pool_(thread_pool_, std::move(params))
    {
        // Run the pool in the background (this is performed by the thread_pool thread).
        // When sync_pool is destroyed, this task will be stopped and joined automatically.
        conn_pool_.async_run(boost::asio::detached);
    }

    // Retrieves a connection from the pool. Throws an exception on error
    boost::mysql::pooled_connection get_connection(
        std::chrono::steady_clock::duration timeout = std::chrono::seconds(30)
    )
    {
        // use_future returns a std::future<pooled_connection>.
        // Calling get() waits for the future to complete and throws an exception on failure.
        // with_diagnostics ensures that the exception contains any server-supplied information.
        // cancel_after applies a timeout to the operation
        // <-
        // clang-format off
        // ->
        return conn_pool_
            .async_get_connection(
                with_diagnostics(boost::asio::cancel_after(timeout, boost::asio::use_future))
            )
            .get();
        // <-
        // clang-format on
        // ->
    }
};
//]

BOOST_AUTO_TEST_CASE(section_connection_pool)
{
    auto server_hostname = get_hostname();

    {
        //[connection_pool_create
        // pool_params contains configuration for the pool.
        // You must specify enough information to establish a connection,
        // including the server address and credentials.
        // You can configure a lot of other things, like pool limits
        boost::mysql::pool_params params;
        params.server_address.emplace_host_and_port(server_hostname);
        params.username = mysql_username;
        params.password = mysql_password;
        params.database = "boost_mysql_examples";

        // The I/O context, required by all I/O operations
        boost::asio::io_context ctx;

        // Construct a pool of connections. The context will be used internally
        // to create the connections and other I/O objects
        boost::mysql::connection_pool pool(ctx, std::move(params));

        // You need to call async_run on the pool before doing anything useful with it.
        // async_run creates connections and keeps them healthy. It must be called
        // only once per pool.
        // The detached completion token means that we don't want to be notified when
        // the operation ends. It's similar to a no-op callback.
        pool.async_run(boost::asio::detached);
        //]

#ifdef BOOST_ASIO_HAS_CO_AWAIT
        run_coro(ctx, [&pool]() -> boost::asio::awaitable<void> {
            co_await get_num_employees(pool);
            pool.cancel();
        });
#endif
    }
    {
        boost::asio::io_context ctx;

        //[connection_pool_configure_size
        boost::mysql::pool_params params;

        // Set the usual params
        params.server_address.emplace_host_and_port(server_hostname);
        params.username = mysql_username;
        params.password = mysql_password;
        params.database = "boost_mysql_examples";

        // Create 10 connections at startup, and allow up to 1000 connections
        params.initial_size = 10;
        params.max_size = 1000;

        boost::mysql::connection_pool pool(ctx, std::move(params));
        //]

#ifdef BOOST_ASIO_HAS_CO_AWAIT
        pool.async_run(boost::asio::detached);
        run_coro(ctx, [&pool]() -> boost::asio::awaitable<void> {
            co_await return_without_reset(pool);
            co_await apply_timeout(pool);
            pool.cancel();
        });
#endif
    }
    {
        //[connection_pool_thread_safe
        // The I/O context, required by all I/O operations
        boost::asio::io_context ctx;

        // The usual pool configuration params
        boost::mysql::pool_params params;
        params.server_address.emplace_host_and_port(server_hostname);
        params.username = mysql_username;
        params.password = mysql_password;
        params.database = "boost_mysql_examples";
        params.thread_safe = true;  // enable thread safety

        // Construct a thread-safe pool
        boost::mysql::connection_pool pool(ctx, std::move(params));

        // We can now pass a reference to pool to other threads,
        // and call async_get_connection concurrently without problem.
        // Individual connections are still not thread-safe.
        //]
    }
    {
        boost::mysql::pool_params params;
        params.server_address.emplace_host_and_port(server_hostname);
        params.username = mysql_username;
        params.password = mysql_password;
        params.database = "boost_mysql_examples";

        sync_pool spool(std::move(params));

        auto conn1 = spool.get_connection();
        BOOST_TEST(conn1.valid());
    }
}

}  // namespace
