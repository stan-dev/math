//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/detached.hpp>
#include <boost/asio/thread_pool.hpp>

#include <cstddef>

#include "tsan_pool_common.hpp"

using boost::mysql::error_code;
namespace mysql = boost::mysql;
namespace asio = boost::asio;
using namespace boost::mysql::test;

namespace {

class task
{
    mysql::connection_pool& pool_;
    coordinator& coord_;
    mysql::results r_;
    mysql::diagnostics diag_;
    mysql::pooled_connection conn_;

public:
    task(mysql::connection_pool& pool, coordinator& coord) : pool_(pool), coord_(coord) {}

    void start_get_connection()
    {
        pool_.async_get_connection(diag_, [this](error_code ec, mysql::pooled_connection c) {
            check_ec(ec, diag_);
            conn_ = std::move(c);
            start_execute();
        });
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
    asio::thread_pool ctx(8);
    mysql::connection_pool pool(ctx, create_pool_params(hostname));
    coordinator coord(pool);
    std::vector<task> tasks;

    // Run the pool
    pool.async_run(asio::detached);

    // Create tasks
    for (std::size_t i = 0; i < num_tasks; ++i)
        tasks.emplace_back(pool, coord);

    // Launch
    for (auto& t : tasks)
        t.start_get_connection();

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
