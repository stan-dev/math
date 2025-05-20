/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/redis/logger.hpp>
#include <boost/asio/detached.hpp>
#include <boost/system/errc.hpp>
#define BOOST_TEST_MODULE conn-run-cancel
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "common.hpp"

#ifdef BOOST_ASIO_HAS_CO_AWAIT
#include <boost/asio/experimental/awaitable_operators.hpp>
#include <boost/asio/as_tuple.hpp>

namespace net = boost::asio;

using boost::redis::operation;
using boost::redis::connection;
using boost::system::error_code;
using net::as_tuple;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore;
using boost::redis::logger;
using namespace std::chrono_literals;

using namespace net::experimental::awaitable_operators;

auto async_cancel_run_with_timer() -> net::awaitable<void>
{
   auto ex = co_await net::this_coro::executor;
   connection conn{ex};

   net::steady_timer st{ex};
   st.expires_after(1s);

   error_code ec1, ec2;
   auto cfg = make_test_config();
   logger l;
   co_await (conn.async_run(cfg, l, redir(ec1)) || st.async_wait(redir(ec2)));

   BOOST_CHECK_EQUAL(ec1, boost::asio::error::operation_aborted);
   BOOST_TEST(!ec2);
}

BOOST_AUTO_TEST_CASE(cancel_run_with_timer)
{
   net::io_context ioc;
   net::co_spawn(ioc.get_executor(), async_cancel_run_with_timer(), net::detached);
   ioc.run();
}

auto
async_check_cancellation_not_missed(int n, std::chrono::milliseconds ms) -> net::awaitable<void>
{
   auto ex = co_await net::this_coro::executor;
   connection conn{ex};

   net::steady_timer timer{ex};

   for (auto i = 0; i < n; ++i) {
      timer.expires_after(ms);
      error_code ec1, ec2;
      auto cfg = make_test_config();
      logger l;
      co_await (conn.async_run(cfg, l, redir(ec1)) || timer.async_wait(redir(ec2)));
      BOOST_CHECK_EQUAL(ec1, boost::asio::error::operation_aborted);
      std::cout << "Counter: " << i << std::endl;
   }
}

// See PR #29
BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_0)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(10, std::chrono::milliseconds{0}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_2)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{2}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_8)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{8}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_16)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{16}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_32)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{32}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_64)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{64}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_128)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{128}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_256)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{256}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_512)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{512}), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(check_implicit_cancel_not_missed_1024)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_check_cancellation_not_missed(20, std::chrono::milliseconds{1024}), net::detached);
   ioc.run();
}

#else
BOOST_AUTO_TEST_CASE(dummy)
{
   BOOST_TEST(true);
}
#endif
