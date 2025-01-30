/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/redis/response.hpp>
#include <boost/system/errc.hpp>
#define BOOST_TEST_MODULE check-health
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <thread>
#include "common.hpp"

namespace net = boost::asio;
namespace redis = boost::redis;
using error_code = boost::system::error_code;
using connection = boost::redis::connection;
using boost::redis::request;
using boost::redis::ignore;
using boost::redis::operation;
using boost::redis::generic_response;
using boost::redis::consume_one;

// TODO: Test cancel(health_check) 

struct push_callback {
   connection* conn1;
   connection* conn2;
   generic_response* resp2;
   request* req1;
   int i = 0;
   boost::asio::coroutine coro{};

   void operator()(error_code ec = {}, std::size_t = 0)
   {
      BOOST_ASIO_CORO_REENTER (coro) for (;;)
      {
         BOOST_ASIO_CORO_YIELD
         conn2->async_receive(*this);
         if (ec) {
            std::clog << "Exiting." << std::endl;
            return;
         }

         BOOST_TEST(resp2->has_value());
         BOOST_TEST(!resp2->value().empty());
         std::clog << "Event> " << resp2->value().front().value << std::endl;
         consume_one(*resp2);

         ++i;

         if (i == 5) {
            std::clog << "Pausing the server" << std::endl;
            // Pause the redis server to test if the health-check exits.
            BOOST_ASIO_CORO_YIELD
            conn1->async_exec(*req1, ignore, *this);
            std::clog << "After pausing> " << ec.message() << std::endl;
            // Don't know in CI we are getting: Got RESP3 simple-error.
            //BOOST_TEST(!ec);
            conn2->cancel(operation::run);
            conn2->cancel(operation::receive);
            conn2->cancel(operation::reconnection);
            return;
         }
      }
   };
};

BOOST_AUTO_TEST_CASE(check_health)
{
   net::io_context ioc;
   connection conn1{ioc};

   request req1;
   req1.push("CLIENT", "PAUSE", "10000", "ALL");

   auto cfg1 = make_test_config();
   cfg1.health_check_id = "conn1";
   cfg1.reconnect_wait_interval = std::chrono::seconds::zero();
   error_code res1;
   conn1.async_run(cfg1, {}, [&](auto ec) {
      std::cout << "async_run 1 completed: " << ec.message() << std::endl;
      res1 = ec;
   });

   //--------------------------------

   // It looks like client pause does not work for clients that are
   // sending MONITOR. I will therefore open a second connection.
   connection conn2{ioc};

   auto cfg2 = make_test_config();
   cfg2.health_check_id = "conn2";
   error_code res2;
   conn2.async_run(cfg2, {}, [&](auto ec){
      std::cout << "async_run 2 completed: " << ec.message() << std::endl;
      res2 = ec;
   });

   request req2;
   req2.push("MONITOR");
   generic_response resp2;
   conn2.set_receive_response(resp2);

   conn2.async_exec(req2, ignore, [](auto ec, auto) {
      std::cout << "async_exec: " << std::endl;
      BOOST_TEST(!ec);
   });

   //--------------------------------
   
   push_callback{&conn1, &conn2, &resp2, &req1}(); // Starts reading pushes.

   ioc.run();

   BOOST_TEST(!!res1);
   BOOST_TEST(!!res2);

   // Waits before exiting otherwise it might cause subsequent tests
   // to fail.
   std::this_thread::sleep_for(std::chrono::seconds{10});
}

