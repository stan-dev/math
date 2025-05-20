//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>
#include <boost/mysql/throw_on_error.hpp>

#include <boost/system/system_error.hpp>
#include <boost/test/unit_test.hpp>

#include "test_common/create_diagnostics.hpp"
#include "test_common/printing.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_throw_on_error)

BOOST_AUTO_TEST_CASE(success)
{
    error_code ec;
    auto diag = create_server_diag("abc");
    BOOST_CHECK_NO_THROW(throw_on_error(ec, diag));
}

BOOST_AUTO_TEST_CASE(server_failure)
{
    error_code ec(boost::mysql::common_server_errc::er_bad_db_error);
    auto diag = create_server_diag("abc");
    BOOST_CHECK_EXCEPTION(
        throw_on_error(ec, diag),
        error_with_diagnostics,
        [&](const error_with_diagnostics& err) {
            BOOST_TEST(err.what() == string_view("er_bad_db_error [mysql.common-server:1049]"));
            BOOST_TEST(err.code() == ec);
            BOOST_TEST(err.get_diagnostics() == diag);
            return true;
        }
    );
}

BOOST_AUTO_TEST_CASE(client_failure)
{
    error_code ec(boost::mysql::client_errc::incomplete_message);
    auto diag = create_client_diag("abc");
    BOOST_CHECK_EXCEPTION(
        throw_on_error(ec, diag),
        error_with_diagnostics,
        [&](const error_with_diagnostics& err) {
            BOOST_TEST(
                err.what() ==
                string_view("abc: An incomplete message was received from the server [mysql.client:1]")
            );
            BOOST_TEST(err.code() == ec);
            BOOST_TEST(err.get_diagnostics() == diag);
            return true;
        }
    );
}

BOOST_AUTO_TEST_CASE(client_failure_no_message)
{
    error_code ec(boost::mysql::client_errc::incomplete_message);
    BOOST_CHECK_EXCEPTION(
        throw_on_error(ec, diagnostics()),
        error_with_diagnostics,
        [&](const error_with_diagnostics& err) {
            BOOST_TEST(
                err.what() ==
                string_view("An incomplete message was received from the server [mysql.client:1]")
            );
            BOOST_TEST(err.code() == ec);
            BOOST_TEST(err.get_diagnostics() == diagnostics());
            return true;
        }
    );
}

BOOST_AUTO_TEST_CASE(no_diagnostics)
{
    error_code ec(boost::mysql::client_errc::incomplete_message);
    BOOST_CHECK_EXCEPTION(throw_on_error(ec), error_with_diagnostics, [&](const error_with_diagnostics& err) {
        BOOST_TEST(
            err.what() == string_view("An incomplete message was received from the server [mysql.client:1]")
        );
        BOOST_TEST(err.code() == ec);
        BOOST_TEST(err.get_diagnostics() == diagnostics());
        return true;
    });
}

BOOST_AUTO_TEST_SUITE_END()
