/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/redis/logger.hpp>
#include <boost/system/errc.hpp>
#define BOOST_TEST_MODULE conn-exec-error
#include <boost/test/included/unit_test.hpp>
#include "common.hpp"
#include <iostream>

namespace net = boost::asio;
namespace redis = boost::redis;
namespace resp3 = redis::resp3;
using error_code = boost::system::error_code;
using connection = boost::redis::connection;
using boost::redis::request;
using boost::redis::response;
using boost::redis::generic_response;
using boost::redis::ignore;
using boost::redis::ignore_t;
using boost::redis::error;
using boost::redis::logger;
using boost::redis::operation;
using namespace std::chrono_literals;

BOOST_AUTO_TEST_CASE(no_ignore_error)
{
   request req;

   // HELLO expects a number, by feeding a string we should get a simple error.
   req.push("HELLO", "not-a-number");

   net::io_context ioc;

   auto conn = std::make_shared<connection>(ioc);

   conn->async_exec(req, ignore, [&](auto ec, auto){
      BOOST_CHECK_EQUAL(ec, error::resp3_simple_error);
      conn->cancel(operation::run);
      conn->cancel(operation::reconnection);
   });

   run(conn);

   ioc.run();
}

BOOST_AUTO_TEST_CASE(has_diagnostic)
{
   request req;

   // HELLO expects a number, by feeding a string we should get a simple error.
   req.push("HELLO", "not-a-number");

   // The second command should be also executed. Notice PING does not
   // require resp3.
   req.push("PING", "Barra do Una");

   net::io_context ioc;

   auto conn = std::make_shared<connection>(ioc);

   response<std::string, std::string> resp;
   conn->async_exec(req, resp, [&](auto ec, auto){
      BOOST_TEST(!ec);

      // HELLO
      BOOST_TEST(std::get<0>(resp).has_error());
      BOOST_CHECK_EQUAL(std::get<0>(resp).error().data_type, resp3::type::simple_error);
      auto const diag = std::get<0>(resp).error().diagnostic;
      BOOST_TEST(!std::empty(diag));
      std::cout << "has_diagnostic: " << diag << std::endl;

      // PING
      BOOST_TEST(std::get<1>(resp).has_value());
      BOOST_CHECK_EQUAL(std::get<1>(resp).value(), "Barra do Una");

      conn->cancel(operation::run);
      conn->cancel(operation::reconnection);
   });

   run(conn);

   ioc.run();
}

BOOST_AUTO_TEST_CASE(resp3_error_in_cmd_pipeline)
{
   request req1;
   req1.push("HELLO", "3");
   req1.push("PING", "req1-msg1");
   req1.push("PING", "req1-msg2", "extra arg"); // Error.
   req1.push("PING", "req1-msg3"); // Should run ok.

   response<ignore_t, std::string, std::string, std::string> resp1;

   request req2;
   req2.push("PING", "req2-msg1");

   response<std::string> resp2;

   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   auto c2 = [&](auto ec, auto)
   {
      BOOST_TEST(!ec);
      BOOST_TEST(std::get<0>(resp2).has_value());
      BOOST_CHECK_EQUAL(std::get<0>(resp2).value(), "req2-msg1");
      conn->cancel(operation::run);
      conn->cancel(operation::reconnection);
   };

   auto c1 = [&](auto ec, auto)
   {
      BOOST_TEST(!ec);
      BOOST_TEST(std::get<2>(resp1).has_error());
      BOOST_CHECK_EQUAL(std::get<2>(resp1).error().data_type, resp3::type::simple_error);
      auto const diag = std::get<2>(resp1).error().diagnostic;
      BOOST_TEST(!std::empty(diag));
      std::cout << "resp3_error_in_cmd_pipeline: " << diag << std::endl;

      BOOST_TEST(std::get<3>(resp1).has_value());
      BOOST_CHECK_EQUAL(std::get<3>(resp1).value(), "req1-msg3");

      conn->async_exec(req2, resp2, c2);
   };

   conn->async_exec(req1, resp1, c1);
   run(conn);

   ioc.run();
}

BOOST_AUTO_TEST_CASE(error_in_transaction)
{
   request req;
   req.push("HELLO", 3);
   req.push("MULTI");
   req.push("PING");
   req.push("PING", "msg2", "error"); // Error.
   req.push("PING");
   req.push("EXEC");
   req.push("PING");

   response<
      ignore_t, // hello
      ignore_t, // multi
      ignore_t, // ping
      ignore_t, // ping
      ignore_t, // ping
      response<std::string, std::string, std::string>, // exec
      std::string // ping
   > resp;


   net::io_context ioc;

   auto conn = std::make_shared<connection>(ioc);

   conn->async_exec(req, resp, [&](auto ec, auto){
      BOOST_TEST(!ec);
      
      BOOST_TEST(std::get<0>(resp).has_value());
      BOOST_TEST(std::get<1>(resp).has_value());
      BOOST_TEST(std::get<2>(resp).has_value());
      BOOST_TEST(std::get<3>(resp).has_value());
      BOOST_TEST(std::get<4>(resp).has_value());
      BOOST_TEST(std::get<5>(resp).has_value());

      // Test errors in the pipeline commands.
      BOOST_TEST(std::get<0>(std::get<5>(resp).value()).has_value());
      BOOST_CHECK_EQUAL(std::get<0>(std::get<5>(resp).value()).value(), "PONG");

      // The ping in the transaction that should be an error.
      BOOST_TEST(std::get<1>(std::get<5>(resp).value()).has_error());
      BOOST_CHECK_EQUAL(std::get<1>(std::get<5>(resp).value()).error().data_type, resp3::type::simple_error);
      auto const diag = std::get<1>(std::get<5>(resp).value()).error().diagnostic;
      BOOST_TEST(!std::empty(diag));

      // The ping thereafter in the transaction should not be an error.
      BOOST_TEST(std::get<2>(std::get<5>(resp).value()).has_value());
      //BOOST_CHECK_EQUAL(std::get<2>(std::get<4>(resp).value()).value(), "PONG");

      // The command right after the pipeline should be successful.
      BOOST_TEST(std::get<6>(resp).has_value());
      BOOST_CHECK_EQUAL(std::get<6>(resp).value(), "PONG");

      conn->cancel(operation::run);
      conn->cancel(operation::reconnection);
   });

   run(conn);

   ioc.run();
}

// This test is important because a SUBSCRIBE command has no response
// on success, but does on error. for example when using a wrong
// syntax, the server will send a simple error response the client is
// not expecting.
//
// Sending the subscribe after the ping command below is just a
// convenience to avoid have it merged in a pipeline making things
// even more complex. For example, without a ping, we might get the
// sequence HELLO + SUBSCRIBE + PING where the hello and ping are
// automatically sent by the implementation. In this case, if the
// subscribe synthax is wrong, redis will send a response, which does
// not exist on success. That response will be interprested as the
// response to the PING command that comes thereafter and won't be
// forwarded to the receive_op, resulting in a difficult to handle
// error.
BOOST_AUTO_TEST_CASE(subscriber_wrong_syntax)
{
   request req1;
   req1.push("PING");

   request req2;
   req2.push("SUBSCRIBE"); // Wrong command synthax.

   net::io_context ioc;
   auto conn = std::make_shared<connection>(ioc);

   auto c2 = [&](auto ec, auto)
   {
      std::cout << "async_exec: subscribe" << std::endl;
      BOOST_TEST(!ec);
   };

   auto c1 = [&](auto ec, auto)
   {
      std::cout << "async_exec: hello" << std::endl;
      BOOST_TEST(!ec);
      conn->async_exec(req2, ignore, c2);
   };

   conn->async_exec(req1, ignore, c1);

   generic_response gresp;
   conn->set_receive_response(gresp);

   auto c3 = [&](auto ec, auto)
   {
      std::cout << "async_receive" << std::endl;
      BOOST_TEST(!ec);
      BOOST_TEST(gresp.has_error());
      BOOST_CHECK_EQUAL(gresp.error().data_type, resp3::type::simple_error);
      BOOST_TEST(!std::empty(gresp.error().diagnostic));
      std::cout << gresp.error().diagnostic << std::endl;
      conn->cancel(operation::run);
      conn->cancel(operation::reconnection);
   };

   conn->async_receive(c3);

   run(conn);

   ioc.run();
}

