/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/detached.hpp>
#include <iostream>

namespace asio = boost::asio;
using boost::redis::connection;
using boost::redis::request;
using boost::redis::response;
using boost::redis::config;

auto main(int argc, char * argv[]) -> int
{
   try {
      config cfg;

      if (argc == 3) {
         cfg.addr.host = argv[1];
         cfg.addr.port = argv[2];
      }

      request req;
      req.push("PING", "Hello world");

      response<std::string> resp;

      asio::io_context ioc;
      connection conn{ioc};

      conn.async_run(cfg, {}, asio::detached);

      conn.async_exec(req, resp, [&](auto ec, auto) {
         if (!ec)
            std::cout << "PING: " << std::get<0>(resp).value() << std::endl;
         conn.cancel();
      });

      ioc.run();

   } catch (std::exception const& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;
   }
}

