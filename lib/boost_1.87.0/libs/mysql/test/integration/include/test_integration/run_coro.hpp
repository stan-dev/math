//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_RUN_CORO_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_RUN_CORO_HPP

#include <boost/asio/awaitable.hpp>
#include <boost/asio/io_context.hpp>

#include <functional>

#include "test_common/source_location.hpp"

namespace boost {
namespace mysql {
namespace test {

#ifdef BOOST_ASIO_HAS_CO_AWAIT
void run_coro(
    asio::io_context& ctx,
    std::function<boost::asio::awaitable<void>(void)> fn,
    source_location loc = BOOST_MYSQL_CURRENT_LOCATION
);
#endif

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
