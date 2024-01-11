/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/redis/logger.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/consign.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/asio/signal_set.hpp>
#include <iostream>

#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace asio = boost::asio;
using namespace std::chrono_literals;
using boost::redis::request;
using boost::redis::generic_response;
using boost::redis::consume_one;
using boost::redis::logger;
using boost::redis::config;
using boost::redis::ignore;
using boost::redis::error;
using boost::system::error_code;
using boost::redis::connection;
using signal_set = asio::deferred_t::as_default_on_t<asio::signal_set>;

/* This example will subscribe and read pushes indefinitely.
 *
 * To test send messages with redis-cli
 *
 *    $ redis-cli -3
 *    127.0.0.1:6379> PUBLISH channel some-message
 *    (integer) 3
 *    127.0.0.1:6379>
 *
 * To test reconnection try, for example, to close all clients currently
 * connected to the Redis instance
 *
 * $ redis-cli
 * > CLIENT kill TYPE pubsub
 */

// Receives server pushes.
auto
receiver(std::shared_ptr<connection> conn) -> asio::awaitable<void>
{
   request req;
   req.push("SUBSCRIBE", "channel");

   generic_response resp;
   conn->set_receive_response(resp);

   // Loop while reconnection is enabled
   while (conn->will_reconnect()) {

      // Reconnect to the channels.
      co_await conn->async_exec(req, ignore, asio::deferred);

      // Loop reading Redis pushs messages.
      for (error_code ec;;) {
         // First tries to read any buffered pushes.
         conn->receive(ec);
         if (ec == error::sync_receive_push_failed) {
            ec = {};
            co_await conn->async_receive(asio::redirect_error(asio::use_awaitable, ec));
         }

         if (ec)
            break; // Connection lost, break so we can reconnect to channels.

         std::cout
            << resp.value().at(1).value
            << " " << resp.value().at(2).value
            << " " << resp.value().at(3).value
            << std::endl;

         consume_one(resp);
      }
   }
}

auto co_main(config cfg) -> asio::awaitable<void>
{
   auto ex = co_await asio::this_coro::executor;
   auto conn = std::make_shared<connection>(ex);
   asio::co_spawn(ex, receiver(conn), asio::detached);
   conn->async_run(cfg, {}, asio::consign(asio::detached, conn));

   signal_set sig_set(ex, SIGINT, SIGTERM);
   co_await sig_set.async_wait();

   conn->cancel();
}

#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
