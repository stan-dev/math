/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/redis/logger.hpp>
#include <boost/system/errc.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/as_tuple.hpp>
#define BOOST_TEST_MODULE conn-push
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "common.hpp"

namespace net = boost::asio;
namespace redis = boost::redis;

using boost::redis::operation;
using connection = boost::redis::connection;
using error_code = boost::system::error_code;
using net::as_tuple;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore;
using boost::redis::ignore_t;
using boost::system::error_code;
using boost::redis::logger;
using namespace std::chrono_literals;

BOOST_AUTO_TEST_CASE(receives_push_waiting_resps)
{
   request req1;
   req1.push("HELLO", 3);
   req1.push("PING", "Message1");

   request req2;
   req2.push("SUBSCRIBE", "channel");

   request req3;
   req3.push("PING", "Message2");
   req3.push("QUIT");

   net::io_context ioc;

   auto conn = std::make_shared<connection>(ioc);

   auto c3 =[](auto ec, auto...)
   {
      std::cout << "c3: " << ec.message() << std::endl;
   };

   auto c2 =[&, conn](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req3, ignore, c3);
   };

   auto c1 =[&, conn](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore, c2);
   };

   conn->async_exec(req1, ignore, c1);

   run(conn, make_test_config(), {});

   bool push_received = false;
   conn->async_receive([&, conn](auto ec, auto){
      std::cout << "async_receive" << std::endl;
      BOOST_TEST(!ec);
      push_received = true;
      conn->cancel();
   });

   ioc.run();

   BOOST_TEST(push_received);
}

BOOST_AUTO_TEST_CASE(push_received1)
{
   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   // Trick: Uses SUBSCRIBE because this command has no response or
   // better said, its response is a server push, which is what we
   // want to test. We send two because we want to test both
   // async_receive and receive.
   request req;
   req.push("SUBSCRIBE", "channel1");
   req.push("SUBSCRIBE", "channel2");

   conn->async_exec(req, ignore, [conn](auto ec, auto){
      std::cout << "async_exec" << std::endl;
      BOOST_TEST(!ec);
   });

   bool push_async_received = false;
   conn->async_receive([&, conn](auto ec, auto){
      std::cout << "(1) async_receive" << std::endl;

      BOOST_TEST(!ec);
      push_async_received = true;

      // Receives the second push synchronously.
      error_code ec2;
      std::size_t res = 0;
      res = conn->receive(ec2);
      BOOST_TEST(!ec2);
      BOOST_TEST(res != std::size_t(0));

      // Tries to receive a third push synchronously.
      ec2 = {};
      res = conn->receive(ec2);
      BOOST_CHECK_EQUAL(ec2, boost::redis::make_error_code(boost::redis::error::sync_receive_push_failed));

      conn->cancel();
   });

   run(conn);
   ioc.run();

   BOOST_TEST(push_async_received);
}

BOOST_AUTO_TEST_CASE(push_filtered_out)
{
   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   request req;
   req.push("HELLO", 3);
   req.push("PING");
   req.push("SUBSCRIBE", "channel");
   req.push("QUIT");

   response<ignore_t, std::string, std::string> resp;
   conn->async_exec(req, resp, [conn](auto ec, auto){
      BOOST_TEST(!ec);
   });

   conn->async_receive([&, conn](auto ec, auto){
      BOOST_TEST(!ec);
      conn->cancel(operation::reconnection);
   });

   run(conn);

   ioc.run();

   BOOST_CHECK_EQUAL(std::get<1>(resp).value(), "PONG");
   BOOST_CHECK_EQUAL(std::get<2>(resp).value(), "OK");
}

#ifdef BOOST_ASIO_HAS_CO_AWAIT
net::awaitable<void>
push_consumer1(std::shared_ptr<connection> conn, bool& push_received)
{
   {
      auto [ec, ev] = co_await conn->async_receive(as_tuple(net::use_awaitable));
      BOOST_TEST(!ec);
   }

   {
      auto [ec, ev] = co_await conn->async_receive(as_tuple(net::use_awaitable));
      BOOST_CHECK_EQUAL(ec, boost::system::errc::errc_t::operation_canceled);
   }

   push_received = true;
}

struct response_error_tag{};
response_error_tag error_tag_obj;

struct response_error_adapter {
   void
   operator()(
      std::size_t, boost::redis::resp3::basic_node<std::string_view> const&, boost::system::error_code& ec)
   {
      ec = boost::redis::error::incompatible_size;
   }

   [[nodiscard]]
   auto get_supported_response_size() const noexcept
      { return static_cast<std::size_t>(-1);}
};

auto boost_redis_adapt(response_error_tag&)
{
   return response_error_adapter{};
}

BOOST_AUTO_TEST_CASE(test_push_adapter)
{
   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   request req;
   req.push("HELLO", 3);
   req.push("PING");
   req.push("SUBSCRIBE", "channel");
   req.push("PING");

   conn->set_receive_response(error_tag_obj);

   conn->async_receive([&, conn](auto ec, auto) {
      BOOST_CHECK_EQUAL(ec, boost::asio::experimental::error::channel_cancelled);
      conn->cancel(operation::reconnection);
   });

   conn->async_exec(req, ignore, [](auto ec, auto){
      BOOST_CHECK_EQUAL(ec, boost::system::errc::errc_t::operation_canceled);
   });

   auto cfg = make_test_config();
   conn->async_run(cfg, {}, [](auto ec){
      BOOST_CHECK_EQUAL(ec, redis::error::incompatible_size);
   });

   ioc.run();

   // TODO: Reset the ioc reconnect and send a quit to ensure
   // reconnection is possible after an error.
}

net::awaitable<void> push_consumer3(std::shared_ptr<connection> conn)
{
   for (;;) {
      co_await conn->async_receive(net::use_awaitable);
   }
}

BOOST_AUTO_TEST_CASE(many_subscribers)
{
   request req0;
   req0.get_config().cancel_on_connection_lost = false;
   req0.push("HELLO", 3);

   request req1;
   req1.get_config().cancel_on_connection_lost = false;
   req1.push("PING", "Message1");

   request req2;
   req2.get_config().cancel_on_connection_lost = false;
   req2.push("SUBSCRIBE", "channel");

   request req3;
   req3.get_config().cancel_on_connection_lost = false;
   req3.push("QUIT");

   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   auto c11 =[&](auto ec, auto...)
   {
      std::cout << "quit sent: " << ec.message() << std::endl;
      conn->cancel(operation::reconnection);
   };
   auto c10 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req3, ignore, c11);
   };
   auto c9 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore, c10);
   };
   auto c8 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req1, ignore,  c9);
   };
   auto c7 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore,  c8);
   };
   auto c6 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore,  c7);
   };
   auto c5 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req1, ignore,  c6);
   };
   auto c4 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore,  c5);
   };
   auto c3 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req1, ignore,  c4);
   };
   auto c2 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore,  c3);
   };
   auto c1 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore,  c2);
   };
   auto c0 =[&](auto ec, auto...)
   {
      BOOST_TEST(!ec);
      conn->async_exec(req1, ignore,  c1);
   };

   conn->async_exec(req0, ignore,  c0);

   run(conn, make_test_config(), {});

   net::co_spawn(ioc.get_executor(), push_consumer3(conn), net::detached);
   ioc.run();
}
#endif
