/* Copyright (c) 2018-2024 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/redis/logger.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#define BOOST_TEST_MODULE conn-quit
#include <boost/test/included/unit_test.hpp>
#include <chrono>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include "common.hpp"

namespace net = boost::asio;
using boost::redis::request;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore;
using boost::redis::logger;
using boost::redis::config;
using boost::redis::operation;
using boost::redis::connection;
using boost::system::error_code;
using namespace std::chrono_literals;

BOOST_AUTO_TEST_CASE(issue_181)
{
   using connection_base = boost::redis::detail::connection_base<net::any_io_executor>;

   auto const level = boost::redis::logger::level::debug;
   net::io_context ioc;
   auto ctx = net::ssl::context{net::ssl::context::tlsv12_client};
   connection_base conn{ioc.get_executor(), std::move(ctx), 1000000};
   net::steady_timer timer{ioc};
   timer.expires_after(std::chrono::seconds{1});

   auto run_cont = [&](auto ec){
      std::cout << "async_run1: " << ec.message() << std::endl;
   };

   auto cfg = make_test_config();
   cfg.health_check_interval = std::chrono::seconds{0};
   cfg.reconnect_wait_interval = std::chrono::seconds{0};
   conn.async_run(cfg, boost::redis::logger{level}, run_cont);
   BOOST_TEST(!conn.run_is_canceled());

   // Uses a timer to wait some time until run has been called.
   auto timer_cont = [&](auto ec){
      std::cout << "timer_cont: " << ec.message() << std::endl;
      BOOST_TEST(!conn.run_is_canceled());
      conn.cancel(operation::run);
      BOOST_TEST(conn.run_is_canceled());
   };

   timer.async_wait(timer_cont);

   ioc.run();
}
