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

// NOTE1: I have observed that if hello and
// blpop are sent toguether, Redis will send the response of hello
// right away, not waiting for blpop.

namespace net = boost::asio;
using error_code = boost::system::error_code;
using namespace net::experimental::awaitable_operators;
using boost::redis::operation;
using boost::redis::error;
using boost::redis::request;
using boost::redis::response;
using boost::redis::generic_response;
using boost::redis::ignore;
using boost::redis::ignore_t;
using boost::redis::logger;
using boost::redis::connection;
using namespace std::chrono_literals;

auto implicit_cancel_of_req_written() -> net::awaitable<void>
{
   auto ex = co_await net::this_coro::executor;
   auto conn = std::make_shared<connection>(ex);

   auto cfg = make_test_config();
   cfg.health_check_interval = std::chrono::seconds::zero();
   run(conn, cfg);

   // See NOTE1.
   request req0;
   req0.push("PING");
   co_await conn->async_exec(req0, ignore, net::use_awaitable);

   // Will be cancelled after it has been written but before the
   // response arrives.
   request req1;
   req1.push("BLPOP", "any", 3);

   net::steady_timer st{ex};
   st.expires_after(std::chrono::seconds{1});

   // Achieves implicit cancellation when the timer fires.
   boost::system::error_code ec1, ec2;
   co_await (
      conn->async_exec(req1, ignore, redir(ec1)) ||
      st.async_wait(redir(ec2))
   );

   conn->cancel();

   // I have observed this produces terminal cancellation so it can't
   // be ignored, an error is expected.
   BOOST_CHECK_EQUAL(ec1, net::error::operation_aborted);
   BOOST_TEST(!ec2);
}

BOOST_AUTO_TEST_CASE(test_ignore_implicit_cancel_of_req_written)
{
   net::io_context ioc;
   net::co_spawn(ioc, implicit_cancel_of_req_written(), net::detached);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(test_cancel_of_req_written_on_run_canceled)
{
   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   request req0;
   req0.push("PING");

   // Sends a request that will be blocked forever, so we can test
   // canceling it while waiting for a response.
   request req1;
   req1.get_config().cancel_on_connection_lost = true;
   req1.get_config().cancel_if_unresponded = true;
   req1.push("BLPOP", "any", 0);

   auto c1 = [&](auto ec, auto)
   {
      BOOST_CHECK_EQUAL(ec, net::error::operation_aborted);
   };

   auto c0 = [&](auto ec, auto)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req1, ignore, c1);
   };

   conn->async_exec(req0, ignore, c0);

   auto cfg = make_test_config();
   cfg.health_check_interval = std::chrono::seconds{5};
   run(conn);

   net::steady_timer st{ioc};
   st.expires_after(std::chrono::seconds{1});
   st.async_wait([&](auto ec){
      BOOST_TEST(!ec);
      conn->cancel(operation::run);
      conn->cancel(operation::reconnection);
   });

   ioc.run();
}

#else
BOOST_AUTO_TEST_CASE(dummy)
{
   BOOST_TEST(true);
}
#endif
