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
#include <boost/mysql/results.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/strand.hpp>
#include <boost/asio/thread_pool.hpp>

#include <algorithm>
#include <chrono>
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
    std::chrono::milliseconds timeout_{1};  // make it likely to get some cancellations

public:
    task(mysql::connection_pool& pool, coordinator& coord)
        : pool_(pool), coord_(coord), strand_(asio::make_strand(pool_.get_executor()))
    {
    }

    void start()
    {
        // get_connection should be run within a strand because cancel_after implies
        // a parallel operation. This is orthogonal to the pool's thread-safety
        // (this is strand is per-task, the pool's is shared between all tasks).
        asio::dispatch(asio::bind_executor(strand_, [this]() { start_get_connection(); }));
    }

    void start_get_connection()
    {
        pool_.async_get_connection(
            diag_,
            asio::cancel_after(
                timeout_,
                asio::bind_executor(
                    strand_,
                    [this](error_code ec, mysql::pooled_connection c) {
                        if (ec == boost::mysql::client_errc::no_connection_available)
                        {
                            // We got a cancellation. Increase timeout so we don't get stuck
                            // and try again
                            timeout_ *= 2;
                            start_get_connection();
                        }
                        else
                        {
                            check_ec(ec, diag_);
                            conn_ = std::move(c);
                            start_execute();
                        }
                    }
                )
            )
        );
    }

    void start_execute()
    {
        auto callback = [this](error_code ec) {
            check_ec(ec, diag_);
            conn_ = boost::mysql::pooled_connection();
            if (coord_.on_iteration_finish())
            {
                start_get_connection();
            }
        };
        conn_->async_execute("SELECT 1", r_, diag_, asio::bind_executor(strand_, std::move(callback)));
    }

    bool had_cancellations() const { return timeout_.count() > 1; }
};

void run(const char* hostname)
{
    // Setup
    asio::thread_pool ctx(8);
    mysql::pool_params params;
    mysql::connection_pool pool(ctx, create_pool_params(hostname));
    coordinator coord(pool);
    std::vector<task> tasks;

    // Run the pool
    pool.async_run(asio::detached);

    // Create connections
    for (std::size_t i = 0; i < num_tasks; ++i)
        tasks.emplace_back(pool, coord);

    // Launch
    for (auto& t : tasks)
        t.start();

    // Run
    ctx.join();

    // Verify that at least one task had a cancellation
    if (!std::any_of(tasks.begin(), tasks.end(), [](const task& t) { return t.had_cancellations(); }))
    {
        std::cerr << "No task had any cancellations\n";
        exit(1);
    }
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
