/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#define BOOST_TEST_MODULE run
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "common.hpp"

namespace net = boost::asio;
namespace redis = boost::redis;

using connection = redis::connection;
using redis::config;
using redis::logger;
using redis::operation;
using boost::system::error_code;
using namespace std::chrono_literals;

bool is_host_not_found(error_code ec)
{
   if (ec == net::error::netdb_errors::host_not_found) return true;
   if (ec == net::error::netdb_errors::host_not_found_try_again) return true;
   return false;
}

BOOST_AUTO_TEST_CASE(resolve_bad_host)
{
   net::io_context ioc;

   auto cfg = make_test_config();
   cfg.addr.host = "Atibaia";
   cfg.addr.port = "6379";
   cfg.resolve_timeout = 10h;
   cfg.connect_timeout = 10h;
   cfg.health_check_interval = 10h;
   cfg.reconnect_wait_interval = 0s;

   auto conn = std::make_shared<connection>(ioc);
   conn->async_run(cfg, {}, [](auto ec){
      BOOST_TEST(is_host_not_found(ec));
   });

   ioc.run();
}

BOOST_AUTO_TEST_CASE(resolve_with_timeout)
{
   net::io_context ioc;

   auto cfg = make_test_config();
   cfg.addr.host = "occase.de";
   cfg.addr.port = "6379";
   cfg.resolve_timeout = 1ms;
   cfg.connect_timeout = 1ms;
   cfg.health_check_interval = 10h;
   cfg.reconnect_wait_interval = 0s;

   auto conn = std::make_shared<connection>(ioc);
   run(conn, cfg);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(connect_bad_port)
{
   net::io_context ioc;

   auto cfg = make_test_config();
   cfg.addr.host = "127.0.0.1";
   cfg.addr.port = "1";
   cfg.resolve_timeout = 10h;
   cfg.connect_timeout = 10s;
   cfg.health_check_interval = 10h;
   cfg.reconnect_wait_interval = 0s;

   auto conn = std::make_shared<connection>(ioc);
   run(conn, cfg, net::error::connection_refused);
   ioc.run();
}

// Hard to test.
//BOOST_AUTO_TEST_CASE(connect_with_timeout)
//{
//   net::io_context ioc;
//
//   config cfg;
//   cfg.addr.host = "example.com";
//   cfg.addr.port = "80";
//   cfg.resolve_timeout = 10s;
//   cfg.connect_timeout = 1ns;
//   cfg.health_check_interval = 10h;
//
//   auto conn = std::make_shared<connection>(ioc);
//   run(conn, cfg, boost::redis::error::connect_timeout);
//   ioc.run();
//}

