/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/consign.hpp>
#include <iostream>

#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace asio = boost::asio;
using boost::redis::request;
using boost::redis::response;
using boost::redis::config;
using boost::redis::connection;

// Called from the main function (see main.cpp)
auto co_main(config cfg) -> asio::awaitable<void>
{
   auto conn = std::make_shared<connection>(co_await asio::this_coro::executor);
   conn->async_run(cfg, {}, asio::consign(asio::detached, conn));

   // A request containing only a ping command.
   request req;
   req.push("PING", "Hello world");

   // Response where the PONG response will be stored.
   response<std::string> resp;

   // Executes the request.
   co_await conn->async_exec(req, resp, asio::deferred);
   conn->cancel();

   std::cout << "PING: " << std::get<0>(resp).value() << std::endl;
}

#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
