//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/socket_stream.hpp>

#include <boost/asio/buffered_stream.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ip/udp.hpp>
#include <boost/asio/local/stream_protocol.hpp>
#include <boost/asio/ssl/stream.hpp>
#include <boost/asio/strand.hpp>
#include <boost/asio/windows/stream_handle.hpp>

using namespace boost::mysql::detail;
namespace asio = boost::asio;

namespace {

struct stream_archetype
{
    using lowest_layer_type = asio::ip::tcp::socket;
    lowest_layer_type& lowest_layer() noexcept
    {
        static asio::io_context ctx;
        static lowest_layer_type res{ctx};
        return res;
    }
};

struct stream_bad_type
{
    using lowest_layer_type = asio::ip::udp::socket;
    lowest_layer_type& lowest_layer() noexcept
    {
        static asio::io_context ctx;
        static lowest_layer_type res{ctx};
        return res;
    }
};

// The streams we regularly use are accepted
static_assert(is_socket_stream<asio::ip::tcp::socket>::value, "");
#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
static_assert(is_socket_stream<asio::local::stream_protocol::socket>::value, "");
#endif
static_assert(is_socket_stream<asio::ssl::stream<asio::ip::tcp::socket>>::value, "");

// Regular streams with more exotic arguments are also accepted
static_assert(is_socket_stream<asio::ssl::stream<asio::ip::tcp::socket&>>::value, "");
static_assert(
    is_socket_stream<asio::ip::tcp::socket::rebind_executor<asio::io_context::executor_type>::other>::value,
    ""
);

// Having several layers works
static_assert(is_socket_stream<asio::buffered_stream<asio::ip::tcp::socket>>::value, "");
static_assert(is_socket_stream<asio::ssl::stream<asio::buffered_stream<asio::ip::tcp::socket>>>::value, "");

// A minimal archetype is accepted
static_assert(is_socket_stream<stream_archetype>::value, "");

// Bad lowest_layer_type
static_assert(!is_socket_stream<stream_bad_type>::value, "");

// Stream that is not a socket
#ifdef BOOST_ASIO_HAS_WINDOWS_STREAM_HANDLE
static_assert(!is_socket_stream<asio::windows::stream_handle>::value, "");
#endif

// No lowest_layer_type
static_assert(!is_socket_stream<int>::value, "");

}  // namespace
