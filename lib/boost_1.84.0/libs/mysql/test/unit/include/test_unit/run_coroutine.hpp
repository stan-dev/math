//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_RUN_COROUTINE_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_RUN_COROUTINE_HPP

#include <boost/asio/use_awaitable.hpp>

#ifdef BOOST_ASIO_HAS_CO_AWAIT

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/use_future.hpp>

#include "test_common/netfun_helpers.hpp"

namespace boost {
namespace mysql {
namespace test {

inline void run_coroutine(asio::any_io_executor ex, std::function<asio::awaitable<void>(void)> coro)
{
    auto fut = boost::asio::co_spawn(ex, std::move(coro), boost::asio::use_future);
    run_until_completion(ex);
    fut.get();
}
}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif  // BOOST_ASIO_HAS_CO_AWAIT

#endif
