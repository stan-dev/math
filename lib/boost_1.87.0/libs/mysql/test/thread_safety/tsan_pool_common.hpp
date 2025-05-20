//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_THREAD_SAFETY_TSAN_POOL_COMMON_HPP
#define BOOST_MYSQL_TEST_THREAD_SAFETY_TSAN_POOL_COMMON_HPP

#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/pool_params.hpp>
#include <boost/mysql/ssl_mode.hpp>

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <mutex>

namespace boost {
namespace mysql {
namespace test {

// Exit on error (thread-safe)
inline void check_ec(error_code ec, diagnostics& diag)
{
    if (ec)
    {
        static std::mutex mtx;
        std::lock_guard<std::mutex> guard{mtx};
        std::cerr << ec << ", " << diag.server_message() << std::endl;
        exit(1);
    }
}

// Constants
constexpr int num_tasks = 100;
constexpr int iterations_per_task = 20;

// Coordinate tasks
class coordinator
{
    std::atomic_int remaining_queries_{num_tasks * iterations_per_task};
    std::atomic_int outstanding_tasks_{num_tasks};
    connection_pool& pool_;

    void on_task_finish()
    {
        if (--outstanding_tasks_ == 0)
            pool_.cancel();
    }

public:
    coordinator(connection_pool& pool) : pool_(pool) {}

    // Call after a task finishes an iteration. Returns true if we should continue
    bool on_iteration_finish()
    {
        bool should_continue = (--remaining_queries_ > 0);
        if (!should_continue)
            on_task_finish();
        return should_continue;
    }
};

// Default pool params
inline pool_params create_pool_params(const char* hostname, std::size_t initial_size = 1)
{
    pool_params params;
    params.server_address.emplace_host_and_port(hostname);
    params.username = "integ_user";
    params.password = "integ_password";
    params.initial_size = initial_size;
    params.max_size = num_tasks - 10;       // create some contention
    params.ssl = mysql::ssl_mode::require;  // double check sharing SSL contexts is OK
    params.thread_safe = true;
    return params;
}

inline void usage(const char* progname)
{
    std::cerr << "Usage: " << progname << " <hostname>\n";
    exit(1);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
