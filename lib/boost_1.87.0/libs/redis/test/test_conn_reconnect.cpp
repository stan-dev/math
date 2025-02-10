/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/detached.hpp>
#define BOOST_TEST_MODULE conn-reconnect
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "common.hpp"

#ifdef BOOST_ASIO_HAS_CO_AWAIT
#include <boost/asio/experimental/awaitable_operators.hpp>

namespace net = boost::asio;
using boost::system::error_code;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore;
using boost::redis::logger;
using boost::redis::operation;
using boost::redis::connection;
using namespace std::chrono_literals;

using namespace boost::asio::experimental::awaitable_operators;

net::awaitable<void> test_reconnect_impl()
{
   auto ex = co_await net::this_coro::executor;

   request req;
   req.push("QUIT");

   auto conn = std::make_shared<connection>(ex);
   run(conn);

   int i = 0;
   for (; i < 5; ++i) {
      error_code ec1, ec2;
      auto cfg = make_test_config();
      logger l;
      co_await conn->async_exec(req, ignore, net::redirect_error(net::use_awaitable, ec1));
      //BOOST_TEST(!ec);
      std::cout << "test_reconnect: " << i << " " << ec2.message() << " " << ec1.message() << std::endl;
   }

   conn->cancel();
   BOOST_CHECK_EQUAL(i, 5);
   co_return;
}

// Test whether the client works after a reconnect.
BOOST_AUTO_TEST_CASE(test_reconnect)
{
   net::io_context ioc;
   net::co_spawn(ioc, test_reconnect_impl(), net::detached);
   ioc.run();
}

auto async_test_reconnect_timeout() -> net::awaitable<void>
{
   auto ex = co_await net::this_coro::executor;

   net::steady_timer st{ex};

   auto conn = std::make_shared<connection>(ex);
   error_code ec1, ec3;

   request req1;
   req1.get_config().cancel_if_not_connected = false;
   req1.get_config().cancel_on_connection_lost = true;
   req1.get_config().cancel_if_unresponded = true;
   req1.push("BLPOP", "any", 0);

   st.expires_after(std::chrono::seconds{1});
   auto cfg = make_test_config();
   co_await (
      conn->async_exec(req1, ignore, redir(ec1)) ||
      st.async_wait(redir(ec3))
   );

   //BOOST_TEST(!ec1);
   //BOOST_TEST(!ec3);

   request req2;
   req2.get_config().cancel_if_not_connected = false;
   req2.get_config().cancel_on_connection_lost = true;
   req2.get_config().cancel_if_unresponded= true;
   req2.push("QUIT");

   st.expires_after(std::chrono::seconds{1});
   co_await (
      conn->async_exec(req1, ignore, net::redirect_error(net::use_awaitable, ec1)) ||
      st.async_wait(net::redirect_error(net::use_awaitable, ec3))
   );
   conn->cancel();

   std::cout << "ccc" << std::endl;

   BOOST_CHECK_EQUAL(ec1, boost::asio::error::operation_aborted);
}

BOOST_AUTO_TEST_CASE(test_reconnect_and_idle)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_test_reconnect_timeout(), net::detached);
   ioc.run();
}
#else
BOOST_AUTO_TEST_CASE(dummy)
{
   BOOST_TEST(true);
}
#endif
