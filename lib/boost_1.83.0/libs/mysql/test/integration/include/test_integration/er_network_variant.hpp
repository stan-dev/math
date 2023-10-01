//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_ER_NETWORK_VARIANT_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_ER_NETWORK_VARIANT_HPP

#include <boost/mysql/string_view.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/core/span.hpp>

#include "test_integration/er_connection.hpp"

namespace boost {
namespace mysql {
namespace test {

class er_network_variant
{
public:
    virtual ~er_network_variant() {}
    virtual bool supports_ssl() const = 0;
    virtual bool is_unix_socket() const = 0;
    virtual const char* stream_name() const = 0;
    virtual const char* variant_name() const = 0;
    std::string name() const { return std::string(stream_name()) + '_' + variant_name(); }
    virtual er_connection_ptr create_connection(boost::asio::any_io_executor, boost::asio::ssl::context&) = 0;
};

boost::span<er_network_variant*> all_variants();
er_network_variant* get_variant(string_view name);

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
