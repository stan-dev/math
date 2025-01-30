//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_ANY_CONNECTION_FIXTURE_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_ANY_CONNECTION_FIXTURE_HPP

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/connect_params.hpp>

#include <boost/assert/source_location.hpp>

#include "test_common/io_context_fixture.hpp"
#include "test_common/source_location.hpp"

namespace boost {
namespace mysql {
namespace test {

struct any_connection_fixture : io_context_fixture
{
    any_connection conn;

    any_connection_fixture() : any_connection_fixture(any_connection_params{}) {}
    any_connection_fixture(any_connection_params params);
    any_connection_fixture(asio::ssl::context& ssl_ctx);
    ~any_connection_fixture();

    void connect(const connect_params& params, source_location loc = BOOST_MYSQL_CURRENT_LOCATION);
    void connect(source_location loc = BOOST_MYSQL_CURRENT_LOCATION);
    void start_transaction(source_location loc = BOOST_MYSQL_CURRENT_LOCATION);
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
