//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_STREAMS_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_STREAMS_HPP

#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/local/stream_protocol.hpp>
#include <boost/asio/ssl/stream.hpp>
#include <boost/asio/use_future.hpp>

namespace boost {
namespace mysql {
namespace test {

// The actual streams we will be using to test
using tcp_socket = boost::asio::ip::tcp::socket;
using tcp_ssl_socket = boost::asio::ssl::stream<tcp_socket>;

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
using unix_socket = boost::asio::local::stream_protocol::socket;
using unix_ssl_socket = boost::asio::ssl::stream<unix_socket>;
#endif

// Stream names
template <class Stream>
constexpr const char* get_stream_name();
template <>
constexpr const char* get_stream_name<tcp_socket>()
{
    return "tcp";
}
template <>
constexpr const char* get_stream_name<tcp_ssl_socket>()
{
    return "tcp_ssl";
}

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
template <>
constexpr const char* get_stream_name<unix_socket>()
{
    return "unix";
}
template <>
constexpr const char* get_stream_name<unix_ssl_socket>()
{
    return "unix_ssl";
}
#endif

// Supports SSL (doesn't use the lib's type trait for test independance)
template <class Stream>
constexpr bool supports_ssl();
template <>
constexpr bool supports_ssl<tcp_socket>()
{
    return false;
}
template <>
constexpr bool supports_ssl<tcp_ssl_socket>()
{
    return true;
}

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
template <>
constexpr bool supports_ssl<unix_socket>()
{
    return false;
}
template <>
constexpr bool supports_ssl<unix_ssl_socket>()
{
    return true;
}
#endif

// is_unix_socket
template <class Stream>
constexpr bool is_unix_socket()
{
    return false;
}

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
template <>
constexpr bool is_unix_socket<unix_socket>()
{
    return true;
}
template <>
constexpr bool is_unix_socket<unix_ssl_socket>()
{
    return true;
}
#endif

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif