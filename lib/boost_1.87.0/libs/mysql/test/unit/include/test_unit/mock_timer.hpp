//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_MOCK_TIMER_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_MOCK_TIMER_HPP

#include <boost/asio/basic_waitable_timer.hpp>
#include <boost/asio/wait_traits.hpp>

#include <chrono>

namespace boost {
namespace mysql {
namespace test {

// Like steady_clock, but the current time can be set within tests
struct mock_clock
{
    using rep = typename std::chrono::steady_clock::rep;
    using period = typename std::chrono::steady_clock::period;
    using duration = typename std::chrono::steady_clock::duration;
    using time_point = typename std::chrono::steady_clock::time_point;
    static constexpr bool is_steady = true;
    static time_point now();
    static void advance_time_by(std::chrono::steady_clock::duration dur);
};

// Helper typedef
using mock_timer = asio::basic_waitable_timer<mysql::test::mock_clock>;

}  // namespace test
}  // namespace mysql
}  // namespace boost

namespace boost {
namespace asio {

// Specialization suitable to use for tests. Instructs Asio to
// create physical timers that wait for a zero duration, effectively
// causing Asio's timer service to poll for ready handlers
template <>
struct wait_traits<mysql::test::mock_clock>
{
    static std::chrono::steady_clock::duration to_wait_duration(std::chrono::steady_clock::duration)
    {
        return std::chrono::nanoseconds(0);
    }

    static std::chrono::steady_clock::duration to_wait_duration(std::chrono::steady_clock::time_point)
    {
        return std::chrono::nanoseconds(0);
    }
};

}  // namespace asio
}  // namespace boost

#endif
