//
// echo_server.cpp
// ~~~~~~~~~~~~~~~
//
// Copyright (c) 2003-2022 Christopher M. Kohlhoff (chris at kohlhoff dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <cstdio>
#include <iostream>
#include <boost/asio.hpp>
#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace net = boost::asio;
namespace this_coro = net::this_coro;
using net::ip::tcp;
using net::detached;
using executor_type = net::io_context::executor_type;
using socket_type = net::basic_stream_socket<net::ip::tcp, executor_type>;
using tcp_socket = net::use_awaitable_t<executor_type>::as_default_on_t<socket_type>;
using acceptor_type = net::basic_socket_acceptor<net::ip::tcp, executor_type>;
using tcp_acceptor = net::use_awaitable_t<executor_type>::as_default_on_t<acceptor_type>;
using awaitable_type = net::awaitable<void, executor_type>;
constexpr net::use_awaitable_t<executor_type> use_awaitable;

awaitable_type echo(tcp_socket socket)
{
  try {
     char data[1024];
     for (;;) {
        std::size_t n = co_await socket.async_read_some(net::buffer(data), use_awaitable);
        co_await async_write(socket, net::buffer(data, n), use_awaitable);
     }
  } catch (std::exception const&) {
     //std::printf("echo Exception: %s\n", e.what());
  }
}

awaitable_type listener()
{
  auto ex = co_await this_coro::executor;
  tcp_acceptor acceptor(ex, {tcp::v4(), 55555});
  for (;;) {
     tcp_socket socket = co_await acceptor.async_accept(use_awaitable);
     co_spawn(ex, echo(std::move(socket)), detached);
  }
}

int main()
{
  try {
     net::io_context io_context{BOOST_ASIO_CONCURRENCY_HINT_UNSAFE_IO};
     co_spawn(io_context, listener(), detached);
     io_context.run();
  } catch (std::exception const& e) {
     std::printf("Exception: %s\n", e.what());
  }
}
#else // defined(BOOST_ASIO_HAS_CO_AWAIT)
auto main() -> int {std::cout << "Requires coroutine support." << std::endl; return 1;}
#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
