/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/asio/co_spawn.hpp>
#include <iostream>

#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace asio = boost::asio;
using tcp_socket = asio::deferred_t::as_default_on_t<asio::ip::tcp::socket>;
using tcp_acceptor = asio::deferred_t::as_default_on_t<asio::ip::tcp::acceptor>;
using signal_set = asio::deferred_t::as_default_on_t<asio::signal_set>;
using boost::redis::request;
using boost::redis::response;
using boost::redis::config;
using boost::system::error_code;
using boost::redis::connection;
using namespace std::chrono_literals;

auto echo_server_session(tcp_socket socket, std::shared_ptr<connection> conn) -> asio::awaitable<void>
{
   request req;
   response<std::string> resp;

   for (std::string buffer;;) {
      auto n = co_await asio::async_read_until(socket, asio::dynamic_buffer(buffer, 1024), "\n");
      req.push("PING", buffer);
      co_await conn->async_exec(req, resp, asio::deferred);
      co_await asio::async_write(socket, asio::buffer(std::get<0>(resp).value()));
      std::get<0>(resp).value().clear();
      req.clear();
      buffer.erase(0, n);
   }
}

// Listens for tcp connections.
auto listener(std::shared_ptr<connection> conn) -> asio::awaitable<void>
{
   try {
      auto ex = co_await asio::this_coro::executor;
      tcp_acceptor acc(ex, {asio::ip::tcp::v4(), 55555});
      for (;;)
         asio::co_spawn(ex, echo_server_session(co_await acc.async_accept(), conn), asio::detached);
   } catch (std::exception const& e) {
      std::clog << "Listener: " << e.what() << std::endl;
   }
}

// Called from the main function (see main.cpp)
auto co_main(config cfg) -> asio::awaitable<void>
{
   auto ex = co_await asio::this_coro::executor;
   auto conn = std::make_shared<connection>(ex);
   asio::co_spawn(ex, listener(conn), asio::detached);
   conn->async_run(cfg, {}, asio::consign(asio::detached, conn));

   signal_set sig_set(ex, SIGINT, SIGTERM);
   co_await sig_set.async_wait();
   conn->cancel();
}

#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
