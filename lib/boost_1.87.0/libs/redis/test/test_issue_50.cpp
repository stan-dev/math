/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

// Must come before any asio header, otherwise build fails on msvc.

#include <boost/redis/connection.hpp>
#include <boost/redis/logger.hpp>
#include <boost/asio/as_tuple.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/consign.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#define BOOST_TEST_MODULE conn-quit
#include <boost/test/included/unit_test.hpp>
#include <tuple>
#include <iostream>
#include "common.hpp"

#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace net = boost::asio;
using steady_timer = net::use_awaitable_t<>::as_default_on_t<net::steady_timer>;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore;
using boost::redis::logger;
using boost::redis::config;
using boost::redis::operation;
using boost::redis::connection;
using boost::system::error_code;
using boost::asio::use_awaitable;
using boost::asio::redirect_error;
using namespace std::chrono_literals;

// Push consumer
auto
receiver(std::shared_ptr<connection> conn) -> net::awaitable<void>
{
   std::cout << "uuu" << std::endl;
   while (conn->will_reconnect()) {
      std::cout << "dddd" << std::endl;
      // Loop reading Redis pushs messages.
      for (;;) {
         std::cout << "aaaa" << std::endl;
         error_code ec;
         co_await conn->async_receive(redirect_error(use_awaitable, ec));
         if (ec) {
            std::cout << "Error in async_receive" << std::endl;
            break;
         }
      }
   }

   std::cout << "Exiting the receiver." << std::endl;
}

auto
periodic_task(std::shared_ptr<connection> conn) -> net::awaitable<void>
{
  net::steady_timer timer{co_await net::this_coro::executor};
  for (int i = 0; i < 10; ++i) {
    std::cout << "In the loop: " << i << std::endl;
    timer.expires_after(std::chrono::milliseconds(50));
    co_await timer.async_wait(net::use_awaitable);

    // Key is not set so it will cause an error since we are passing
    // an adapter that does not accept null, this will cause an error
    // that result in the connection being closed.
    request req;
    req.push("GET", "mykey");
    auto [ec, u] = co_await conn->async_exec(req, ignore, net::as_tuple(net::use_awaitable));
    if (ec) {
      std::cout << "(1)Error: " << ec << std::endl;
    } else {
      std::cout << "no error: " << std::endl;
    }
  }

  std::cout << "Periodic task done!" << std::endl;
  conn->cancel(operation::run);
  conn->cancel(operation::receive);
  conn->cancel(operation::reconnection);
}

auto co_main(config) -> net::awaitable<void>
{
   auto ex = co_await net::this_coro::executor;
   auto conn = std::make_shared<connection>(ex);

   net::co_spawn(ex, receiver(conn), net::detached);
   net::co_spawn(ex, periodic_task(conn), net::detached);
   auto cfg = make_test_config();
   conn->async_run(cfg, {}, net::consign(net::detached, conn));
}

BOOST_AUTO_TEST_CASE(issue_50)
{
   net::io_context ioc;
   net::co_spawn(ioc, std::move(co_main({})), net::detached);
   ioc.run();
}

#else // defined(BOOST_ASIO_HAS_CO_AWAIT)

BOOST_AUTO_TEST_CASE(issue_50)
{
}

#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
