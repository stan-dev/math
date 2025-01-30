//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/connection.hpp>
#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/coroutine.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/io_context.hpp>

#include <chrono>
#include <cstddef>
#include <iostream>

using boost::mysql::error_code;
using std::chrono::steady_clock;
namespace mysql = boost::mysql;
namespace asio = boost::asio;

namespace {

static constexpr std::size_t num_parallel = 100;
static constexpr std::size_t total = num_parallel * 100;
static constexpr const char* default_unix_path = "/var/run/mysqld/mysqld.sock";

class coordinator
{
    bool finished_{};
    std::size_t remaining_queries_{total};
    std::size_t outstanding_tasks_{num_parallel};
    steady_clock::time_point tp_start_;
    steady_clock::time_point tp_finish_;
    mysql::connection_pool* pool_{};

public:
    coordinator(mysql::connection_pool* pool = nullptr) : pool_(pool) {}
    std::chrono::milliseconds ellapsed() const
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(tp_finish_ - tp_start_);
    }
    void record_start() { tp_start_ = steady_clock::now(); }
    void on_finish()
    {
        if (--outstanding_tasks_ == 0)
        {
            tp_finish_ = steady_clock::now();
            if (pool_)
                pool_->cancel();
        }
    }
    bool on_loop_finish()
    {
        if (--remaining_queries_ == 0)
            finished_ = true;
        return !finished_;
    }

    bool check_ec(error_code ec, const mysql::diagnostics& diag)
    {
        if (ec)
        {
            finished_ = true;
            std::cout << ec << ", " << diag.server_message() << std::endl;
        }
        return !finished_;
    }
};

class task_nopool
{
    mysql::any_connection conn_;
    mysql::results r_;
    mysql::diagnostics diag_;
    coordinator* coord_{};
    const mysql::connect_params* params_;
    mysql::statement stmt_;
    asio::coroutine coro_;

public:
    task_nopool(asio::any_io_executor ex, coordinator& coord, const mysql::connect_params& params)
        : conn_(ex), coord_(&coord), params_(&params)
    {
    }

    void resume(error_code ec = {})
    {
        // Error checking
        if (!coord_->check_ec(ec, diag_))
        {
            coord_->on_finish();
            return;
        }

        BOOST_ASIO_CORO_REENTER(coro_)
        {
            while (true)
            {
                BOOST_ASIO_CORO_YIELD
                conn_.async_connect(*params_, diag_, [this](error_code ec) { resume(ec); });

                BOOST_ASIO_CORO_YIELD
                conn_.async_prepare_statement(
                    "SELECT tax_id FROM company WHERE id = ?",
                    diag_,
                    [this](error_code ec, boost::mysql::statement s) {
                        stmt_ = s;
                        resume(ec);
                    }
                );

                BOOST_ASIO_CORO_YIELD
                conn_.async_execute(stmt_.bind("HGS"), r_, diag_, [this](error_code ec) { resume(ec); });

                BOOST_ASIO_CORO_YIELD
                conn_.async_close(diag_, [this](error_code ec) { resume(ec); });

                if (!coord_->on_loop_finish())
                {
                    coord_->on_finish();
                    return;
                }
            }
        }
    }
};

class task_pool
{
    mysql::connection_pool* pool_;
    mysql::results r_;
    mysql::diagnostics diag_;
    coordinator* coord_{};
    mysql::pooled_connection conn_;
    asio::coroutine coro_;
    mysql::statement stmt_;

    void on_finish()
    {
        conn_ = mysql::pooled_connection();
        coord_->on_finish();
    }

public:
    task_pool(mysql::connection_pool& pool, coordinator& coord) : pool_(&pool), coord_(&coord) {}

    void resume(error_code ec = {})
    {
        // Error checking
        if (!coord_->check_ec(ec, diag_))
        {
            on_finish();
            return;
        }

        BOOST_ASIO_CORO_REENTER(coro_)
        {
            while (true)
            {
                BOOST_ASIO_CORO_YIELD
                pool_->async_get_connection(diag_, [this](error_code ec, mysql::pooled_connection c) {
                    conn_ = std::move(c);
                    resume(ec);
                });

                BOOST_ASIO_CORO_YIELD
                conn_->async_prepare_statement(
                    "SELECT tax_id FROM company WHERE id = ?",
                    diag_,
                    [this](error_code ec, boost::mysql::statement s) {
                        stmt_ = s;
                        resume(ec);
                    }
                );

                BOOST_ASIO_CORO_YIELD
                conn_->async_execute(stmt_.bind("HGS"), r_, diag_, [this](error_code ec) { resume(ec); });

                conn_ = boost::mysql::pooled_connection();

                if (!coord_->on_loop_finish())
                {
                    on_finish();
                    return;
                }
            }
        }
    }
};

void run_nopool(mysql::any_address server_addr, bool use_ssl)
{
    // Setup
    asio::io_context ctx;
    mysql::connect_params params;
    params.server_address = std::move(server_addr);
    params.username = "example_user";
    params.password = "example_password";
    params.database = "boost_mysql_examples";
    params.ssl = use_ssl ? mysql::ssl_mode::require : mysql::ssl_mode::disable;
    std::vector<task_nopool> conns;
    coordinator coord;

    // Create connections
    for (std::size_t i = 0; i < num_parallel; ++i)
        conns.emplace_back(ctx.get_executor(), coord, params);

    // Launch
    coord.record_start();
    for (auto& conn : conns)
        conn.resume(error_code());

    // Run
    ctx.run();

    // Print ellapsed time
    std::cout << coord.ellapsed().count() << std::flush;
}

void run_pool(mysql::any_address server_addr, bool use_ssl)
{
    // Setup
    asio::io_context ctx;
    mysql::pool_params params;
    params.server_address = std::move(server_addr);
    params.username = "example_user";
    params.password = "example_password";
    params.database = "boost_mysql_examples";
    params.max_size = num_parallel;
    params.ssl = use_ssl ? mysql::ssl_mode::require : mysql::ssl_mode::disable;

    mysql::connection_pool pool(ctx, std::move(params));
    pool.async_run(asio::detached);

    std::vector<task_pool> conns;
    coordinator coord(&pool);

    // Create connections
    for (std::size_t i = 0; i < num_parallel; ++i)
        conns.emplace_back(pool, coord);

    // Launch
    coord.record_start();
    for (auto& conn : conns)
        conn.resume(error_code());

    // Run
    ctx.run();

    // Print ellapsed time
    std::cout << coord.ellapsed().count() << std::flush;
}

static constexpr const char* options[] = {
    "nopool-tcp",
    "nopool-tcpssl",
    "nopool-unix",
    "pool-tcp",
    "pool-tcpssl",
    "pool-unix",
};

void usage(const char* progname)
{
    std::cerr << "Usage: " << progname << " <benchmark-type> <server-addr>\nAvailable options:\n";
    for (const char* opt : options)
        std::cerr << "    " << opt << "\n";
    exit(1);
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        usage(argv[0]);
    }

    mysql::string_view opt = argv[1];
    mysql::string_view addr = argv[2];
    boost::mysql::host_and_port tcp_addr;

    if (opt == "nopool-tcp")
    {
        tcp_addr.host = addr;
        run_nopool(std::move(tcp_addr), false);
    }
    else if (opt == "nopool-tcpssl")
    {
        tcp_addr.host = addr;
        run_nopool(std::move(tcp_addr), true);
    }
    else if (opt == "nopool-unix")
    {
        run_nopool(mysql::unix_path{default_unix_path}, false);
    }
    else if (opt == "pool-tcp")
    {
        tcp_addr.host = addr;
        run_pool(std::move(tcp_addr), false);
    }
    else if (opt == "pool-tcpssl")
    {
        tcp_addr.host = addr;
        run_pool(std::move(tcp_addr), true);
    }
    else if (opt == "pool-unix")
    {
        run_pool(mysql::unix_path{default_unix_path}, false);
    }
    else
        usage(argv[0]);
}
