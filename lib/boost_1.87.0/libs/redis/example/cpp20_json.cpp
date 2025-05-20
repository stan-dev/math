/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/connection.hpp>
#include <boost/asio/deferred.hpp>
#include <boost/asio/detached.hpp>
#include <boost/describe.hpp>
#include <boost/asio/consign.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <string>
#include <iostream>

#if defined(BOOST_ASIO_HAS_CO_AWAIT)

#include <boost/json/serialize.hpp>
#include <boost/json/parse.hpp>
#include <boost/json/value_from.hpp>
#include <boost/json/value_to.hpp>
#include <boost/redis/resp3/serialization.hpp>

namespace asio = boost::asio;
using namespace boost::describe;
using boost::redis::request;
using boost::redis::response;
using boost::redis::ignore_t;
using boost::redis::config;
using boost::redis::connection;

// Struct that will be stored in Redis using json serialization. 
struct user {
   std::string name;
   std::string age;
   std::string country;
};

// The type must be described for serialization to work.
BOOST_DESCRIBE_STRUCT(user, (), (name, age, country))

// Boost.Redis customization points (example/json.hpp)
void boost_redis_to_bulk(std::string& to, user const& u)
   { boost::redis::resp3::boost_redis_to_bulk(to, boost::json::serialize(boost::json::value_from(u))); }

void boost_redis_from_bulk(user& u, std::string_view sv, boost::system::error_code&)
   { u = boost::json::value_to<user>(boost::json::parse(sv)); }

auto co_main(config cfg) -> asio::awaitable<void>
{
   auto ex = co_await asio::this_coro::executor;
   auto conn = std::make_shared<connection>(ex);
   conn->async_run(cfg, {}, asio::consign(asio::detached, conn));

   // user object that will be stored in Redis in json format.
   user const u{"Joao", "58", "Brazil"};

   // Stores and retrieves in the same request.
   request req;
   req.push("SET", "json-key", u); // Stores in Redis.
   req.push("GET", "json-key"); // Retrieves from Redis.

   response<ignore_t, user> resp;

   co_await conn->async_exec(req, resp, asio::deferred);
   conn->cancel();

   // Prints the first ping
   std::cout
      << "Name: " << std::get<1>(resp).value().name << "\n"
      << "Age: " << std::get<1>(resp).value().age << "\n"
      << "Country: " << std::get<1>(resp).value().country << "\n";
}

#endif // defined(BOOST_ASIO_HAS_CO_AWAIT)
