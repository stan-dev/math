//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_categories.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/mariadb_server_errc.hpp>
#include <boost/mysql/mysql_server_errc.hpp>
#include <boost/mysql/results.hpp>

#include "test_integration/tcp_network_fixture.hpp"

using namespace boost::mysql::test;
using boost::mysql::diagnostics;
using boost::mysql::error_code;
using boost::mysql::results;

namespace {

BOOST_AUTO_TEST_SUITE(test_db_specific)

BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mysql5"))
BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mariadb"))
BOOST_FIXTURE_TEST_CASE(mysql_specific_error_code, tcp_network_fixture)
{
    connect();
    error_code ec;
    diagnostics diag;
    results result;

    // This is reported as a common, less desriptive error in MySQL5 and MariaDB
    conn.execute("select * from one_row_table where field_varchar regexp '(('", result, ec, diag);
    error_code expected_ec(
        boost::mysql::mysql_server_errc::er_regexp_mismatched_paren,
        boost::mysql::get_mysql_server_category()
    );
    BOOST_TEST(ec == expected_ec);
}

BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mysql5"))
BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mysql8"))
BOOST_FIXTURE_TEST_CASE(mariadb_specific_error_code, tcp_network_fixture)
{
    connect();
    error_code ec;
    diagnostics diag;
    results result;

    // This is reported as a common error in MySQL5 and MySQL8
    conn.execute("WITH abc AS (SELECT 1), abc as (SELECT 2) SELECT * FROM abc", result, ec, diag);
    error_code expected_ec(
        boost::mysql::mariadb_server_errc::er_dup_query_name,
        boost::mysql::get_mariadb_server_category()
    );
    BOOST_TEST(ec == expected_ec);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
