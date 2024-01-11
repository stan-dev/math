/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <iostream>

#define BOOST_TEST_MODULE request
#include <boost/test/included/unit_test.hpp>

#include <boost/redis/request.hpp>

using boost::redis::request;

// TODO: Serialization.

BOOST_AUTO_TEST_CASE(single_arg_allocator)
{
   request req1;
   req1.push("PING");
   BOOST_CHECK_EQUAL(req1.payload(), std::string{"*1\r\n$4\r\nPING\r\n"});
}

BOOST_AUTO_TEST_CASE(arg_int)
{
   request req;
   req.push("PING", 42);
   BOOST_CHECK_EQUAL(req.payload(), std::string{"*2\r\n$4\r\nPING\r\n$2\r\n42\r\n"});
}

BOOST_AUTO_TEST_CASE(multiple_args)
{
   char const* res = "*5\r\n$3\r\nSET\r\n$3\r\nkey\r\n$5\r\nvalue\r\n$2\r\nEX\r\n$1\r\n2\r\n";
   request req;
   req.push("SET", "key", "value", "EX", "2");
   BOOST_CHECK_EQUAL(req.payload(), std::string{res});
}

BOOST_AUTO_TEST_CASE(container_and_range)
{
   std::map<std::string, std::string> in{{"key1", "value1"}, {"key2", "value2"}};

   char const* res = "*6\r\n$4\r\nHSET\r\n$3\r\nkey\r\n$4\r\nkey1\r\n$6\r\nvalue1\r\n$4\r\nkey2\r\n$6\r\nvalue2\r\n";

   request req1;
   req1.push_range("HSET", "key", in);
   BOOST_CHECK_EQUAL(req1.payload(), std::string{res});

   request req2;
   req2.push_range("HSET", "key", std::cbegin(in), std::cend(in));
   BOOST_CHECK_EQUAL(req2.payload(), std::string{res});
}
