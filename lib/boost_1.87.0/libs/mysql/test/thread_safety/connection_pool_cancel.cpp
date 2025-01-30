//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/pool_params.hpp>

#include <boost/asio/bind_executor.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

#include <memory>

#include "tsan_pool_common.hpp"

using boost::mysql::error_code;
namespace mysql = boost::mysql;
namespace asio = boost::asio;
using namespace boost::mysql::test;

namespace {

void run(const char* hostname)
{
    // Setup
    asio::thread_pool ctx(8);

    for (int i = 0; i < 20; ++i)
    {
        // Create a pool
        auto pool = std::make_shared<mysql::connection_pool>(ctx, create_pool_params(hostname, 10));

        // Run the pool within the thread pool
        asio::post(asio::bind_executor(ctx.get_executor(), [pool]() { pool->async_run(asio::detached); }));

        // Issue a cancellation
        asio::post(asio::bind_executor(ctx.get_executor(), [pool]() { pool->cancel(); }));
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
