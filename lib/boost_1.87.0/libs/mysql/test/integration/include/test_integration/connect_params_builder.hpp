//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_CONNECT_PARAMS_BUILDER_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_CONNECT_PARAMS_BUILDER_HPP

#include <boost/mysql/any_address.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/string_view.hpp>

#include <cstdint>

#include "test_common/ci_server.hpp"

namespace boost {
namespace mysql {
namespace test {

class connect_params_builder
{
    handshake_params res_{integ_user, integ_passwd, integ_db};
    any_address addr_;

public:
    connect_params_builder() { addr_.emplace_host_and_port(get_hostname()); }

    connect_params_builder& set_unix()
    {
        addr_.emplace_unix_path(default_unix_path);
        return *this;
    }

    connect_params_builder& credentials(string_view username, string_view passwd)
    {
        res_.set_username(username);
        res_.set_password(passwd);
        return *this;
    }

    connect_params_builder& database(string_view db)
    {
        res_.set_database(db);
        return *this;
    }

    connect_params_builder& disable_ssl() { return ssl(ssl_mode::disable); }

    connect_params_builder& ssl(ssl_mode ssl)
    {
        res_.set_ssl(ssl);
        return *this;
    }

    connect_params_builder& multi_queries(bool v)
    {
        res_.set_multi_queries(v);
        return *this;
    }

    connect_params_builder& collation(std::uint16_t v)
    {
        res_.set_connection_collation(v);
        return *this;
    }

    handshake_params build_hparams() const { return res_; }

    connect_params build();
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
