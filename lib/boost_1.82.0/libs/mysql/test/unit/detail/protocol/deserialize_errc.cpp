//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/protocol/deserialize_errc.hpp>

#include <boost/test/unit_test.hpp>

using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::detail::deserialize_errc;
using boost::mysql::detail::to_error_code;

namespace {

BOOST_AUTO_TEST_SUITE(test_deserialize_errc)

BOOST_AUTO_TEST_CASE(to_error_code_)
{
    BOOST_TEST(to_error_code(deserialize_errc::ok) == error_code());
    BOOST_TEST(
        to_error_code(deserialize_errc::incomplete_message) ==
        error_code(client_errc::incomplete_message)
    );
    BOOST_TEST(
        to_error_code(deserialize_errc::protocol_value_error) ==
        error_code(client_errc::protocol_value_error)
    );
    BOOST_TEST(
        to_error_code(deserialize_errc::server_unsupported) ==
        error_code(client_errc::server_unsupported)
    );
}

BOOST_AUTO_TEST_SUITE_END()  // test_deserialize_errc

}  // namespace
