//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_TEST_ANY_CONNECTION_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_TEST_ANY_CONNECTION_HPP

#include <boost/mysql/any_connection.hpp>

#include <boost/asio/io_context.hpp>

#include "test_unit/test_stream.hpp"

namespace boost {
namespace mysql {
namespace test {

any_connection create_test_any_connection(asio::io_context& ctx, any_connection_params params = {});
test_stream& get_stream(any_connection& conn);

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
