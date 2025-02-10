/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <iostream>
#include <boost/asio.hpp>
#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace net = boost::asio;

using net::ip::tcp;
using tcp_socket = net::use_awaitable_t<>::as_default_on_t<net::ip::tcp::socket>;
using timer_type = net::use_awaitable_t<>::as_default_on_t<net::steady_timer>;

auto example(boost::asio::ip::tcp::endpoint ep, std::string msg, int n) -> net::awaitable<void>
{
   try {
      auto ex = co_await net::this_coro::executor;

      tcp_socket socket{ex};
      co_await socket.async_connect(ep);

      std::string buffer;
      auto dbuffer = net::dynamic_buffer(buffer);
      for (int i = 0; i < n; ++i) {
         co_await net::async_write(socket, net::buffer(msg));
         auto n = co_await net::async_read_until(socket, dbuffer, "\n");
         //std::printf("> %s", buffer.data());
         dbuffer.consume(n);
      }

      //std::printf("Ok: %s", msg.data());
   } catch (std::exception const& e) {
      std::cerr << "Error: " << e.what() << std::endl;
   }
}

int main(int argc, char* argv[])
{
   try {
      int sessions = 1;
      int msgs = 1;

      if (argc == 3) {
         sessions = std::stoi(argv[1]);
         msgs = std::stoi(argv[2]);
      }

      net::io_context ioc;

      tcp::resolver resv{ioc};
      auto const res = resv.resolve("127.0.0.1", "55555");
      auto ep = *std::begin(res);

      for (int i = 0; i < sessions; ++i)
         net::co_spawn(ioc, example(ep, "Some message\n", msgs), net::detached);

      ioc.run();
   } catch (std::exception const& e) {
      std::cerr << e.what() << std::endl;
   }
}
#else // defined(BOOST_ASIO_HAS_CO_AWAIT)
auto main() -> int {std::cout << "Requires coroutine support." << std::endl; return 1;}
#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
