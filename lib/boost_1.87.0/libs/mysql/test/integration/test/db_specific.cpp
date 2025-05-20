//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/error_categories.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/mariadb_server_errc.hpp>
#include <boost/mysql/mysql_server_errc.hpp>
#include <boost/mysql/results.hpp>

#include "test_common/create_diagnostics.hpp"
#include "test_common/network_result.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/server_features.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_SUITE(test_db_specific)

BOOST_TEST_DECORATOR(*run_if(&server_features::regex_error_codes))
BOOST_FIXTURE_TEST_CASE(mysql_specific_error_code, any_connection_fixture)
{
    connect();
    results result;

    // This is reported as a common, less descriptive error in MySQL5 and MariaDB
    error_code expected_ec(mysql_server_errc::er_regexp_mismatched_paren, get_mysql_server_category());
    conn.async_execute("select * from one_row_table where field_varchar regexp '(('", result, as_netresult)
        .validate_error(expected_ec, create_server_diag("Mismatched parenthesis in regular expression."));
}

BOOST_TEST_DECORATOR(*run_if(&server_features::dup_query_error_codes))
BOOST_FIXTURE_TEST_CASE(mariadb_specific_error_code, any_connection_fixture)
{
    connect();
    results result;

    // This is reported as a common error in MySQL5 and MySQL8
    error_code expected_ec(mariadb_server_errc::er_dup_query_name, get_mariadb_server_category());
    conn.async_execute("WITH abc AS (SELECT 1), abc as (SELECT 2) SELECT * FROM abc", result, as_netresult)
        .validate_error(expected_ec, create_server_diag("Duplicate query name `abc` in WITH clause"));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
