//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_SNIPPETS_FIXTURE_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_SNIPPETS_SNIPPETS_FIXTURE_HPP

#include <boost/mysql/connect_params.hpp>

#include "test_common/ci_server.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/snippets/credentials.hpp"

namespace boost {
namespace mysql {
namespace test {

inline connect_params snippets_connect_params()
{
    connect_params params;
    params.server_address.emplace_host_and_port(get_hostname());
    params.username = mysql_username;
    params.password = mysql_password;
    params.database = "boost_mysql_examples";
    params.ssl = ssl_mode::disable;
    params.multi_queries = true;
    return params;
}

struct snippets_fixture : any_connection_fixture
{
    snippets_fixture() { connect(snippets_connect_params()); }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
