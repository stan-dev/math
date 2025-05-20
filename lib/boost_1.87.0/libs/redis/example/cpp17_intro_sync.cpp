/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include "sync_connection.hpp"

#include <string>
#include <iostream>

using boost::redis::sync_connection;
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

      sync_connection conn;
      conn.run(cfg);

      request req;
      req.push("PING");

      response<std::string> resp;

      conn.exec(req, resp);
      conn.stop();

      std::cout << "Response: " << std::get<0>(resp).value() << std::endl;

   } catch (std::exception const& e) {
      std::cerr << e.what() << std::endl;
   }
}
