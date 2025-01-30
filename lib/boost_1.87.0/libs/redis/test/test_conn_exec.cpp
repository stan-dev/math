/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/system/errc.hpp>
#include <boost/asio/detached.hpp>
#define BOOST_TEST_MODULE conn-exec
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "common.hpp"

// TODO: Test whether HELLO won't be inserted passt commands that have
// been already writen.
// TODO: Test async_exec with empty request e.g. hgetall with an empty
// container.

namespace net = boost::asio;
using boost::redis::config;
using boost::redis::connection;
using boost::redis::generic_response;
using boost::redis::ignore;
using boost::redis::operation;
using boost::redis::request;
using boost::redis::response;

// Sends three requests where one of them has a hello with a priority
// set, which means it should be executed first.
BOOST_AUTO_TEST_CASE(hello_priority)
{
   request req1;
   req1.push("PING", "req1");

   request req2;
   req2.get_config().hello_with_priority = false;
   req2.push("HELLO", 3);
   req2.push("PING", "req2");

   request req3;
   req3.get_config().hello_with_priority = true;
   req3.push("HELLO", 3);
   req3.push("PING", "req3");

   net::io_context ioc;

   auto conn = std::make_shared<connection>(ioc);

   bool seen1 = false;
   bool seen2 = false;
   bool seen3 = false;

   conn->async_exec(req1, ignore, [&](auto ec, auto){
      // Second callback to the called.
      std::cout << "req1" << std::endl;
      BOOST_TEST(!ec);
      BOOST_TEST(!seen2);
      BOOST_TEST(seen3);
      seen1 = true;
   });

   conn->async_exec(req2, ignore, [&](auto ec, auto){
      // Last callback to the called.
      std::cout << "req2" << std::endl;
      BOOST_TEST(!ec);
      BOOST_TEST(seen1);
      BOOST_TEST(seen3);
      seen2 = true;
      conn->cancel(operation::run);
      conn->cancel(operation::reconnection);
   });

   conn->async_exec(req3, ignore, [&](auto ec, auto){
      // Callback that will be called first.
      std::cout << "req3" << std::endl;
      BOOST_TEST(!ec);
      BOOST_TEST(!seen1);
      BOOST_TEST(!seen2);
      seen3 = true;
   });

   run(conn);
   ioc.run();
}

// Tries to receive a string in an int and gets an error.
BOOST_AUTO_TEST_CASE(wrong_response_data_type)
{
   request req;
   req.push("PING");

   // Wrong data type.
   response<int> resp;
   net::io_context ioc;

   auto conn = std::make_shared<connection>(ioc);

   conn->async_exec(req, resp, [conn](auto ec, auto){
      BOOST_CHECK_EQUAL(ec, boost::redis::error::not_a_number);
      conn->cancel(operation::reconnection);
   });

   run(conn);
   ioc.run();
}

BOOST_AUTO_TEST_CASE(cancel_request_if_not_connected)
{
   request req;
   req.get_config().cancel_if_not_connected = true;
   req.push("PING");

   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);
   conn->async_exec(req, ignore, [conn](auto ec, auto){
      BOOST_CHECK_EQUAL(ec, boost::redis::error::not_connected);
      conn->cancel();
   });

   ioc.run();
}

BOOST_AUTO_TEST_CASE(correct_database)
{
   auto cfg = make_test_config();
   cfg.database_index = 2;

   net::io_context ioc;

   auto conn = std::make_shared<connection>(ioc);

   request req;
   req.push("CLIENT", "LIST");

   generic_response resp;

   conn->async_exec(req, resp, [&](auto ec, auto n){
         BOOST_TEST(!ec);
         std::clog << "async_exec has completed: " << n << std::endl;
         conn->cancel();
   });

   conn->async_run(cfg, {}, [](auto){
         std::clog << "async_run has exited." << std::endl;
   });

   ioc.run();

   assert(!resp.value().empty());
   auto const& value = resp.value().front().value;
   auto const pos = value.find("db=");
   auto const index_str = value.substr(pos + 3, 1);
   auto const index = std::stoi(index_str);

   // This check might fail if more than one client is connected to
   // redis when the CLIENT LIST command is run. 
   BOOST_CHECK_EQUAL(cfg.database_index.value(), index);
}

BOOST_AUTO_TEST_CASE(large_number_of_concurrent_requests_issue_170)
{
   // See https://github.com/boostorg/redis/issues/170

   std::string payload;
   payload.resize(1024);
   std::fill(std::begin(payload), std::end(payload), 'A');

   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   auto cfg = make_test_config();
   cfg.health_check_interval = std::chrono::seconds(0);
   conn->async_run(cfg, {}, net::detached);

   int counter = 0;
   int const repeat = 8000;

   for (int i = 0; i < repeat; ++i) {
      auto req = std::make_shared<request>();
      req->push("PING", payload);
      conn->async_exec(*req, ignore, [req, &counter, conn](auto ec, auto) {
         BOOST_TEST(!ec);
         if (++counter == repeat)
            conn->cancel();
      });
   }

   ioc.run();

   BOOST_CHECK_EQUAL(counter, repeat);
}

