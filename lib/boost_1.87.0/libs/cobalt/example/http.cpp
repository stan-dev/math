// Copyright (c) 2023 Matthijs Möhlmann
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/asio.hpp>
#include <boost/asio/ssl.hpp>
#include <boost/asio/ssl/stream_base.hpp>
#include <boost/asio/system_timer.hpp>
#include <boost/cobalt.hpp>
#include <boost/cobalt/promise.hpp>
#include <boost/cobalt/race.hpp>
#include <boost/cobalt/this_thread.hpp>
#include <boost/beast.hpp>

#include <boost/beast/core/flat_buffer.hpp>
#include <boost/beast/http/empty_body.hpp>
#include <boost/beast/http/message.hpp>
#include <boost/beast/http/read.hpp>
#include <boost/beast/http/string_body.hpp>
#include <boost/beast/http/verb.hpp>
#include <boost/beast/websocket/stream.hpp>
#include <stdexcept>

namespace cobalt = boost::cobalt;
namespace beast = boost::beast;

using executor_type = cobalt::use_op_t::executor_with_default<cobalt::executor>;
using socket_type = typename boost::asio::ip::tcp::socket::rebind_executor<
    executor_type>::other;
using ssl_socket_type = boost::asio::ssl::stream<socket_type>;
using acceptor_type = typename boost::asio::ip::tcp::acceptor::rebind_executor<
    executor_type>::other;
using websocket_type = beast::websocket::stream<ssl_socket_type>;


cobalt::promise<ssl_socket_type> connect(std::string_view host,
                                        boost::asio::ssl::context &ctx) {
  boost::asio::ip::tcp::resolver resolve{cobalt::this_thread::get_executor()};
  auto endpoints = co_await resolve.async_resolve(host, "https", cobalt::use_op);

  // Timer for timeouts

  ssl_socket_type sock{cobalt::this_thread::get_executor(), ctx};
  printf("connecting\n");

  co_await sock.next_layer().async_connect(*endpoints.begin());
  printf("connected\n");

  // Connected, now do the handshake
  printf("handshaking\n");
  co_await sock.async_handshake(boost::asio::ssl::stream_base::client);
  printf("hand shook\n");
  co_return sock;
}

cobalt::main co_main(int argc, char **argv)
{
  boost::asio::ssl::context ctx{boost::asio::ssl::context::tls_client};
  auto conn = co_await connect("boost.org", ctx);
  printf("connected\n");
  beast::http::request<beast::http::empty_body> req{beast::http::verb::get, "/index.html", 11};
  req.set(beast::http::field::host, "boost.org");
  co_await beast::http::async_write(conn, req, cobalt::use_op);

  // read the response
  beast::flat_buffer b;
  beast::http::response<beast::http::string_body> response;
  co_await beast::http::async_read(conn, b, response, cobalt::use_op);

  // write the response
  printf("%s\n", response.body().c_str());
  co_return 0;
}
