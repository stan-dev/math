/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/consign.hpp>
#include <iostream>

#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace asio = boost::asio;
using boost::redis::request;
using boost::redis::response;
using boost::redis::config;
using boost::redis::logger;
using boost::redis::connection;

auto verify_certificate(bool, asio::ssl::verify_context&) -> bool
{
   std::cout << "set_verify_callback" << std::endl;
   return true;
}

auto co_main(config cfg) -> asio::awaitable<void>
{
   cfg.use_ssl = true;
   cfg.username = "aedis";
   cfg.password = "aedis";
   cfg.addr.host = "db.occase.de";
   cfg.addr.port = "6380";

   auto conn = std::make_shared<connection>(co_await asio::this_coro::executor);
   conn->async_run(cfg, {}, asio::consign(asio::detached, conn));

   request req;
   req.push("PING");

   response<std::string> resp;

   conn->next_layer().set_verify_mode(asio::ssl::verify_peer);
   conn->next_layer().set_verify_callback(verify_certificate);

   co_await conn->async_exec(req, resp, asio::deferred);
   conn->cancel();

   std::cout << "Response: " << std::get<0>(resp).value() << std::endl;
}

#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
