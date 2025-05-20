/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/system/errc.hpp>
#define BOOST_TEST_MODULE echo-stress
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "common.hpp"

#ifdef BOOST_ASIO_HAS_CO_AWAIT

namespace net = boost::asio;
using error_code = boost::system::error_code;
using boost::redis::operation;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore;
using boost::redis::ignore_t;
using boost::redis::logger;
using boost::redis::connection;
using boost::redis::usage;
using boost::redis::error;

std::ostream& operator<<(std::ostream& os, usage const& u)
{
   os
      << "Commands sent: " << u.commands_sent << "\n"
      << "Bytes sent: " << u.bytes_sent << "\n"
      << "Responses received: " << u.responses_received << "\n"
      << "Pushes received: " << u.pushes_received << "\n"
      << "Response bytes received: " << u.response_bytes_received << "\n"
      << "Push bytes received: " << u.push_bytes_received;

   return os;
}

auto push_consumer(std::shared_ptr<connection> conn, int expected) -> net::awaitable<void>
{
   int c = 0;
   for (error_code ec;;) {
      conn->receive(ec);
      if (ec == error::sync_receive_push_failed) {
         ec = {};
         co_await conn->async_receive(redirect_error(net::use_awaitable, ec));
      } else if (!ec) {
         //std::cout << "Skipping suspension." << std::endl;
      }

      if (ec) {
         BOOST_TEST(false);
         std::cout << "push_consumer error: " << ec.message() << std::endl;
         co_return;
      }
      if (++c == expected)
         break;
   }

   conn->cancel();
}

auto
echo_session(
   std::shared_ptr<connection> conn,
   std::shared_ptr<request> pubs,
   int n) -> net::awaitable<void>
{
   for (auto i = 0; i < n; ++i)
      co_await conn->async_exec(*pubs, ignore, net::deferred);
}

auto async_echo_stress(std::shared_ptr<connection> conn) -> net::awaitable<void>
{
   auto ex = co_await net::this_coro::executor;
   auto cfg = make_test_config();
   cfg.health_check_interval = std::chrono::seconds::zero();
   run(conn, cfg,
       boost::asio::error::operation_aborted,
       boost::redis::operation::receive,
       boost::redis::logger::level::crit);

   request req;
   req.push("SUBSCRIBE", "channel");
   co_await conn->async_exec(req, ignore, net::deferred);

   // Number of coroutines that will send pings sharing the same
   // connection to redis.
   int const sessions = 150;

   // The number of pings that will be sent by each session.
   int const msgs = 200;

   // The number of publishes that will be sent by each session with
   // each message.
   int const n_pubs = 25;

   // This is the total number of pushes we will receive.
   int total_pushes = sessions * msgs * n_pubs + 1;

   auto pubs = std::make_shared<request>();
   pubs->push("PING");
   for (int i = 0; i < n_pubs; ++i)
      pubs->push("PUBLISH", "channel", "payload");

   // Op that will consume the pushes counting down until all expected
   // pushes have been received.
   net::co_spawn(ex, push_consumer(conn, total_pushes), net::detached);

   for (int i = 0; i < sessions; ++i) 
      net::co_spawn(ex, echo_session(conn, pubs, msgs), net::detached);
}

BOOST_AUTO_TEST_CASE(echo_stress)
{
   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);
   net::co_spawn(ioc, async_echo_stress(conn), net::detached);
   ioc.run();

   std::cout
      << "-------------------\n"
      << conn->get_usage()
      << std::endl;
}

#else
BOOST_AUTO_TEST_CASE(dummy)
{
   BOOST_TEST(true);
}
#endif
