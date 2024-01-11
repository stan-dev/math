/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/co_spawn.hpp>
#include <map>
#include <vector>
#include <iostream>

#if defined(BOOST_ASIO_HAS_CO_AWAIT)

namespace asio = boost::asio;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore_t;
using boost::redis::ignore;
using boost::redis::config;
using boost::redis::connection;
using boost::asio::awaitable;
using boost::asio::deferred;
using boost::asio::detached;
using boost::asio::consign;

void print(std::map<std::string, std::string> const& cont)
{
   for (auto const& e: cont)
      std::cout << e.first << ": " << e.second << "\n";
}

void print(std::vector<int> const& cont)
{
   for (auto const& e: cont) std::cout << e << " ";
   std::cout << "\n";
}

// Stores the content of some STL containers in Redis.
auto store(std::shared_ptr<connection> conn) -> awaitable<void>
{
   std::vector<int> vec
      {1, 2, 3, 4, 5, 6};

   std::map<std::string, std::string> map
      {{"key1", "value1"}, {"key2", "value2"}, {"key3", "value3"}};

   request req;
   req.push_range("RPUSH", "rpush-key", vec);
   req.push_range("HSET", "hset-key", map);

   co_await conn->async_exec(req, ignore, deferred);
}

auto hgetall(std::shared_ptr<connection> conn) -> awaitable<void>
{
   // A request contains multiple commands.
   request req;
   req.push("HGETALL", "hset-key");

   // Responses as tuple elements.
   response<std::map<std::string, std::string>> resp;

   // Executes the request and reads the response.
   co_await conn->async_exec(req, resp, deferred);

   print(std::get<0>(resp).value());
}

// Retrieves in a transaction.
auto transaction(std::shared_ptr<connection> conn) -> awaitable<void>
{
   request req;
   req.push("MULTI");
   req.push("LRANGE", "rpush-key", 0, -1); // Retrieves
   req.push("HGETALL", "hset-key"); // Retrieves
   req.push("EXEC");

   response<
      ignore_t, // multi
      ignore_t, // lrange
      ignore_t, // hgetall
      response<std::optional<std::vector<int>>, std::optional<std::map<std::string, std::string>>> // exec
   > resp;

   co_await conn->async_exec(req, resp, deferred);

   print(std::get<0>(std::get<3>(resp).value()).value().value());
   print(std::get<1>(std::get<3>(resp).value()).value().value());
}

// Called from the main function (see main.cpp)
awaitable<void> co_main(config cfg)
{
   auto conn = std::make_shared<connection>(co_await asio::this_coro::executor);
   conn->async_run(cfg, {}, consign(detached, conn));

   co_await store(conn);
   co_await transaction(conn);
   co_await hgetall(conn);
   conn->cancel();
}

#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
