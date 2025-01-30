//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_POLL_UNTIL_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_POLL_UNTIL_HPP

#include <boost/asio/io_context.hpp>

#include <functional>

#include "test_common/source_location.hpp"

namespace boost {
namespace mysql {
namespace test {

// poll until *done == true
void poll_until(asio::io_context& ctx, const bool* done, source_location loc = BOOST_MYSQL_CURRENT_LOCATION);

// poll until done() == true
void poll_until(
    asio::io_context& ctx,
    const std::function<bool()>& done,
    source_location loc = BOOST_MYSQL_CURRENT_LOCATION
);

// run fn() in ctx, then poll until it completes
void run_in_context(
    asio::io_context& ctx,
    const std::function<void()>& fn,
    source_location loc = BOOST_MYSQL_CURRENT_LOCATION
);

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
