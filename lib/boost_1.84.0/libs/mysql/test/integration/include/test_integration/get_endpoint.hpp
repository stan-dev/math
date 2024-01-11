//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_GET_ENDPOINT_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_GET_ENDPOINT_HPP

#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/local/stream_protocol.hpp>

#include <stdexcept>
#include <utility>

namespace boost {
namespace mysql {
namespace test {

template <class Protocol>
struct endpoint_getter;

template <>
struct endpoint_getter<boost::asio::ip::tcp>
{
    boost::asio::ip::tcp::endpoint operator()();
};

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
template <>
struct endpoint_getter<boost::asio::local::stream_protocol>
{
    boost::asio::local::stream_protocol::endpoint operator()();
};
#endif

template <class Stream>
typename Stream::lowest_layer_type::endpoint_type get_endpoint()
{
    return endpoint_getter<typename Stream::lowest_layer_type::protocol_type>()();
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
