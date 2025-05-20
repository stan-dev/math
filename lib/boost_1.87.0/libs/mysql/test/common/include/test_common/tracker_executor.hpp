//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_TRACKER_EXECUTOR_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_TRACKER_EXECUTOR_HPP

#include <boost/asio/any_io_executor.hpp>
#include <boost/core/span.hpp>

namespace boost {
namespace mysql {
namespace test {

// A tracker executor wraps an any_io_executor and tracks execute calls.
// Every executor has a distinct ID. When an executor starts executing
// a function, its ID is pushed to a thread-local call stack, and popped
// when the function returns.
// This allows us to reliable check whether "we're running in the context of executor X"
struct tracker_executor_result
{
    int executor_id;
    asio::any_io_executor ex;
};

// Create
tracker_executor_result create_tracker_executor(asio::any_io_executor inner);

// Get the executor call stack, as a span of IDs. Most recent call last.
boost::span<const int> executor_stack();

// Get the ID of a tracker executor, or -1 if it's not a tracker executor
int get_executor_id(asio::any_io_executor);

// We maintain a thread-local state variable that tracks whether we're running
// in the context of an initiation function or not.
// Use this guard when invoking initiation functions to set/clear the flag,
// and the function below to check it.
struct initiation_guard
{
    initiation_guard();
    initiation_guard(const initiation_guard&) = delete;
    initiation_guard(initiation_guard&&) = delete;
    initiation_guard& operator=(const initiation_guard&) = delete;
    initiation_guard& operator=(initiation_guard&&) = delete;
    ~initiation_guard();
};
bool is_initiation_function();

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
