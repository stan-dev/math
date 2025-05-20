//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_LOG_ERROR_HPP
#define BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_LOG_ERROR_HPP

//[example_connection_pool_log_error_hpp
//
// File: log_error.hpp
//

#include <iostream>
#include <mutex>

// Helper function to safely write diagnostics to std::cerr.
// Since we're in a multi-threaded environment, directly writing to std::cerr
// can lead to interleaved output, so we should synchronize calls with a mutex.
// This function is only called in rare cases (e.g. unhandled exceptions),
// so we can afford the synchronization overhead.

namespace notes {

// If you're in C++17+, you can write this using fold expressions
// instead of recursion.
inline void log_error_impl() {}

template <class Arg1, class... Tail>
void log_error_impl(const Arg1& arg, const Tail&... tail)
{
    std::cerr << arg;
    log_error_impl(tail...);
}

template <class... Args>
void log_error(const Args&... args)
{
    static std::mutex mtx;

    // Acquire the mutex, then write the passed arguments to std::cerr.
    std::unique_lock<std::mutex> lock(mtx);
    log_error_impl(args...);
    std::cerr << std::endl;
}

}  // namespace notes

//]

#endif
