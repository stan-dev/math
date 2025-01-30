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

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/thread_pool.hpp>

#include <cstddef>

#include "tsan_pool_common.hpp"

using boost::mysql::error_code;
using namespace boost::mysql::test;
namespace mysql = boost::mysql;
namespace asio = boost::asio;

namespace {

class task
{
    mysql::connection_pool& pool_;
    coordinator& coord_;
    asio::any_io_executor token_ex_;
    mysql::results r_;
    mysql::diagnostics diag_;
    mysql::pooled_connection conn_;

public:
    task(mysql::connection_pool& pool, coordinator& coord, asio::any_io_executor token_ex)
        : pool_(pool), coord_(coord), token_ex_(std::move(token_ex))
    {
    }

    void start_get_connection()
    {
        // Verify that we achieve thread-safety even if the token
        // has a bound executor associated to an execution context different
        // from the one that the pool is using (regression check).
        pool_.async_get_connection(
            diag_,
            asio::bind_executor(
                token_ex_,
                [this](error_code ec, mysql::pooled_connection c) {
                    check_ec(ec, diag_);
                    conn_ = std::move(c);
                    start_execute();
                }
            )
        );
    }

    void start_execute()
    {
        conn_->async_execute("SELECT 1", r_, diag_, [this](error_code ec) {
            check_ec(ec, diag_);
            conn_ = boost::mysql::pooled_connection();
            if (coord_.on_iteration_finish())
            {
                start_get_connection();
            }
        });
    }
};

void run(const char* hostname)
{
    // Setup
    asio::thread_pool pool_ctx(8);
    asio::thread_pool external_ctx(4);
    mysql::connection_pool pool(pool_ctx, create_pool_params(hostname));
    coordinator coord(pool);
    std::vector<task> tasks;

    // The pool should be thread-safe even if we pass a token with a custom
    // executor to async_run
    pool.async_run(asio::bind_executor(external_ctx.get_executor(), asio::detached));

    // Create tasks
    for (std::size_t i = 0; i < num_tasks; ++i)
        tasks.emplace_back(pool, coord, external_ctx.get_executor());

    // Launch
    for (auto& t : tasks)
        t.start_get_connection();

    // Run
    pool_ctx.join();
    external_ctx.join();
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
