/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/system/errc.hpp>
#define BOOST_TEST_MODULE conn-exec-cancel
#include <boost/test/included/unit_test.hpp>
#include <boost/asio/detached.hpp>
#include "common.hpp"
#include <iostream>

#ifdef BOOST_ASIO_HAS_CO_AWAIT
#include <boost/asio/experimental/awaitable_operators.hpp>

// NOTE1: Sends hello separately. I have observed that if hello and
// blpop are sent toguether, Redis will send the response of hello
// right away, not waiting for blpop. That is why we have to send it
// separately.

namespace net = boost::asio;
using error_code = boost::system::error_code;
using namespace net::experimental::awaitable_operators;
using boost::redis::operation;
using boost::redis::request;
using boost::redis::response;
using boost::redis::generic_response;
using boost::redis::ignore;
using boost::redis::ignore_t;
using boost::redis::config;
using boost::redis::logger;
using boost::redis::connection;
using namespace std::chrono_literals;

auto async_ignore_explicit_cancel_of_req_written() -> net::awaitable<void>
{
   auto ex = co_await net::this_coro::executor;

   generic_response gresp;
   auto conn = std::make_shared<connection>(ex);

   run(conn);

   net::steady_timer st{ex};
   st.expires_after(std::chrono::seconds{1});

   // See NOTE1.
   request req0;
   req0.push("PING", "async_ignore_explicit_cancel_of_req_written");
   co_await conn->async_exec(req0, gresp, net::use_awaitable);

   request req1;
   req1.push("BLPOP", "any", 3);

   bool seen = false;
   conn->async_exec(req1, gresp, [&](auto ec, auto) mutable{
      // No error should occur since the cancelation should be
      // ignored.
      std::cout << "async_exec (1): " << ec.message() << std::endl;
      BOOST_TEST(!ec);
      seen = true;
   });

   // Will complete while BLPOP is pending.
   boost::system::error_code ec1;
   co_await st.async_wait(net::redirect_error(net::use_awaitable, ec1));
   conn->cancel(operation::exec);

   BOOST_TEST(!ec1);

   request req3;
   req3.push("PING");

   // Test whether the connection remains usable after a call to
   // cancel(exec).
   co_await conn->async_exec(req3, gresp, net::redirect_error(net::use_awaitable, ec1));
   conn->cancel();

   BOOST_TEST(!ec1);
   BOOST_TEST(seen);
}

BOOST_AUTO_TEST_CASE(test_ignore_explicit_cancel_of_req_written)
{
   net::io_context ioc;
   net::co_spawn(ioc, async_ignore_explicit_cancel_of_req_written(), net::detached);
   ioc.run();
}

#else
BOOST_AUTO_TEST_CASE(dummy)
{
   BOOST_TEST(true);
}
#endif
