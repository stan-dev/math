/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/asio/posix/stream_descriptor.hpp>
#include <unistd.h>
#include <iostream>

#if defined(BOOST_ASIO_HAS_CO_AWAIT)
#if defined(BOOST_ASIO_HAS_POSIX_STREAM_DESCRIPTOR)

namespace asio = boost::asio;
using stream_descriptor = asio::deferred_t::as_default_on_t<asio::posix::stream_descriptor>;
using signal_set = asio::deferred_t::as_default_on_t<asio::signal_set>;
using boost::asio::async_read_until;
using boost::asio::awaitable;
using boost::asio::co_spawn;
using boost::asio::consign;
using boost::asio::deferred;
using boost::asio::detached;
using boost::asio::dynamic_buffer;
using boost::asio::redirect_error;
using boost::asio::use_awaitable;
using boost::redis::config;
using boost::redis::connection;
using boost::redis::generic_response;
using boost::redis::ignore;
using boost::redis::request;
using boost::system::error_code;
using namespace std::chrono_literals;

// Chat over Redis pubsub. To test, run this program from multiple
// terminals and type messages to stdin.

auto
receiver(std::shared_ptr<connection> conn) -> awaitable<void>
{
   request req;
   req.push("SUBSCRIBE", "channel");

   generic_response resp;
   conn->set_receive_response(resp);

   while (conn->will_reconnect()) {

      // Subscribe to channels.
      co_await conn->async_exec(req, ignore, deferred);

      // Loop reading Redis push messages.
      for (error_code ec;;) {
         co_await conn->async_receive(redirect_error(use_awaitable, ec));
         if (ec)
            break; // Connection lost, break so we can reconnect to channels.
         std::cout
            << resp.value().at(1).value
            << " " << resp.value().at(2).value
            << " " << resp.value().at(3).value
            << std::endl;
         resp.value().clear();
      }
   }
}

// Publishes stdin messages to a Redis channel.
auto publisher(std::shared_ptr<stream_descriptor> in, std::shared_ptr<connection> conn) -> awaitable<void>
{
   for (std::string msg;;) {
      auto n = co_await async_read_until(*in, dynamic_buffer(msg, 1024), "\n");
      request req;
      req.push("PUBLISH", "channel", msg);
      co_await conn->async_exec(req, ignore, deferred);
      msg.erase(0, n);
   }
}

// Called from the main function (see main.cpp)
auto co_main(config cfg) -> awaitable<void>
{
   auto ex = co_await asio::this_coro::executor;
   auto conn = std::make_shared<connection>(ex);
   auto stream = std::make_shared<stream_descriptor>(ex, ::dup(STDIN_FILENO));

   co_spawn(ex, receiver(conn), detached);
   co_spawn(ex, publisher(stream, conn), detached);
   conn->async_run(cfg, {}, consign(detached, conn));

   signal_set sig_set{ex, SIGINT, SIGTERM};
   co_await sig_set.async_wait();
   conn->cancel();
   stream->cancel();
}

#else // defined(BOOST_ASIO_HAS_POSIX_STREAM_DESCRIPTOR)
auto co_main(config const&) -> awaitable<void>
{
   std::cout << "Requires support for posix streams." << std::endl;
   co_return;
}
#endif // defined(BOOST_ASIO_HAS_POSIX_STREAM_DESCRIPTOR)
#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
