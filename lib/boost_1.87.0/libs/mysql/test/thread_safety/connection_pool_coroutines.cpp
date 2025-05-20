//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/this_coro.hpp>
#include <boost/asio/thread_pool.hpp>

#include <cstddef>
#include <exception>

#include "tsan_pool_common.hpp"

#ifdef BOOST_ASIO_HAS_CO_AWAIT

namespace mysql = boost::mysql;
namespace asio = boost::asio;
using namespace boost::mysql::test;

namespace {

asio::awaitable<void> task(mysql::connection_pool& pool, coordinator& coord)
{
    mysql::results r;

    do
    {
        // Coroutine handlers always have a bound executor.
        // Verify that we achieve thread-safety even in this case (regression check)
        // Exercise return_without_reset, too
        auto conn = co_await pool.async_get_connection();
        co_await conn->async_execute("SELECT 1", r);
        conn.return_without_reset();
    } while (coord.on_iteration_finish());
}

void rethrow_on_err(std::exception_ptr err)
{
    if (err)
        std::rethrow_exception(err);
}

void run(const char* hostname)
{
    // Setup
    asio::thread_pool ctx(8);
    mysql::pool_params params;
    mysql::connection_pool pool(ctx, create_pool_params(hostname));
    coordinator coord(pool);

    // The pool should be thread-safe even if we pass a token with a custom
    // executor to async_run (as happens with coroutines)
    asio::co_spawn(ctx, [&]() -> asio::awaitable<void> { co_await pool.async_run(); }, rethrow_on_err);

    // Create and launch tasks
    for (std::size_t i = 0; i < num_tasks; ++i)
    {
        asio::co_spawn(ctx, [&]() -> asio::awaitable<void> { return task(pool, coord); }, rethrow_on_err);
    }

    // Run
    ctx.join();
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        usage(argv[0]);
    }

    run(argv[1]);
}

#else

int main() { std::cout << "Sorry, your compiler doesn't support C++20 coroutines\n"; }

#endif
