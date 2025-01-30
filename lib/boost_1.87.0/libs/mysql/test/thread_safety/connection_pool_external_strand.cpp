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

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/strand.hpp>
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
    asio::strand<asio::any_io_executor> strand_;
    mysql::results r_;
    mysql::diagnostics diag_;
    mysql::pooled_connection conn_;

public:
    task(mysql::connection_pool& pool, coordinator& coord, asio::strand<asio::any_io_executor> strand)
        : pool_(pool), coord_(coord), strand_(std::move(strand))
    {
    }

    void start() { enter_strand(); }

    void enter_strand()
    {
        asio::dispatch(asio::bind_executor(strand_, [this]() { start_get_connection(); }));
    }

    void start_get_connection()
    {
        pool_.async_get_connection(diag_, [this](error_code ec, mysql::pooled_connection c) {
            check_ec(ec, diag_);
            conn_ = std::move(c);
            exit_strand();
        });
    }

    void exit_strand()
    {
        asio::post(asio::bind_executor(strand_.get_inner_executor(), [this]() { start_execute(); }));
    }

    void start_execute()
    {
        conn_->async_execute("SELECT 1", r_, diag_, [this](error_code ec) {
            check_ec(ec, diag_);
            asio::dispatch(asio::bind_executor(strand_, [this]() { return_connection(); }));
        });
    }

    void return_connection()
    {
        // Returning a connection must happen within the strand
        conn_ = boost::mysql::pooled_connection();
        if (coord_.on_iteration_finish())
        {
            start_get_connection();
        }
    }
};

// Tests that we can build a thread-safe connection pool ourselves,
// by passing a strand as the pool executor
void run(const char* hostname)
{
    // Setup
    asio::thread_pool ctx(8);
    auto strand = asio::make_strand(ctx);
    auto params = create_pool_params(hostname);
    params.thread_safe = false;
    params.connection_executor = ctx.get_executor();
    mysql::connection_pool pool(strand, std::move(params));
    coordinator coord(pool);
    std::vector<task> tasks;

    // Run the pool
    pool.async_run(asio::detached);

    // Create tasks
    for (std::size_t i = 0; i < num_tasks; ++i)
        tasks.emplace_back(pool, coord, strand);

    // Launch
    for (auto& t : tasks)
        t.start();

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
