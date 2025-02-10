//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/pipeline.hpp>

#include <boost/test/unit_test.hpp>

#include <vector>

#include "test_common/create_basic.hpp"
#include "test_common/create_diagnostics.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
namespace asio = boost::asio;
using boost::test_tools::per_element;

namespace {

BOOST_AUTO_TEST_SUITE(test_pipeline)

BOOST_FIXTURE_TEST_CASE(success, any_connection_fixture)
{
    // Setup
    connect();
    pipeline_request req;
    std::vector<stage_response> res;

    // Populate the request
    req.add_reset_connection()
        .add_execute("SET @myvar = 15")
        .add_prepare_statement("SELECT * FROM one_row_table WHERE id = ?");

    // Run it
    conn.async_run_pipeline(req, res, as_netresult).validate_no_error();

    // Check results
    BOOST_TEST_REQUIRE(res.size() == 3u);
    BOOST_TEST(res[0].error() == error_code());
    BOOST_TEST(res[0].diag() == diagnostics());
    BOOST_TEST(!res[0].has_results());
    BOOST_TEST(!res[0].has_statement());
    BOOST_TEST(res[1].has_results());
    auto stmt = res[2].as_statement();
    BOOST_TEST(stmt.valid());
    BOOST_TEST(conn.current_character_set().error() == client_errc::unknown_character_set);

    // Re-populate the pipeline, with different stages
    req.clear();
    req.add_set_character_set(utf8mb4_charset)
        .add_execute(stmt, 0)
        .add_execute_range(stmt, std::vector<field_view>{field_view(1)})
        .add_close_statement(stmt);

    // Run it
    conn.async_run_pipeline(req, res, as_netresult).validate_no_error();

    // Check results
    BOOST_TEST_REQUIRE(res.size() == 4u);
    BOOST_TEST(res[0].error() == error_code());
    BOOST_TEST(res[0].diag() == diagnostics());
    BOOST_TEST(!res[0].has_results());
    BOOST_TEST(!res[0].has_statement());
    BOOST_TEST(res[1].as_results().rows() == rows(), per_element());
    BOOST_TEST(res[2].as_results().rows() == makerows(2, 1, "f0"), per_element());
    BOOST_TEST(res[3].error() == error_code());
    BOOST_TEST(res[3].diag() == diagnostics());
    BOOST_TEST(!res[3].has_results());
    BOOST_TEST(!res[3].has_statement());
}

BOOST_FIXTURE_TEST_CASE(errors, any_connection_fixture)
{
    // Setup
    connect();
    pipeline_request req;
    std::vector<stage_response> res;

    // Populate the request with some successes and some errors
    req.add_execute("SET @myvar = 42")                                                   // OK
        .add_prepare_statement("SELECT * FROM bad_table WHERE id = ?")                   // error: bad table
        .add_execute("")                                                                 // error: empty query
        .add_execute("SELECT @myvar")                                                    // OK
        .add_set_character_set(character_set{"bad_charset", utf8mb4_charset.next_char})  // bad charset
        .add_execute("SELECT 'abc'");                                                    // OK

    // Run it. The result of the operation is the first encountered error
    conn.async_run_pipeline(req, res, as_netresult)
        .validate_error(
            common_server_errc::er_no_such_table,
            "Table 'boost_mysql_integtests.bad_table' doesn't exist"
        );

    // Check results
    BOOST_TEST_REQUIRE(res.size() == 6u);
    BOOST_TEST(res[0].has_results());
    BOOST_TEST(res[1].error() == common_server_errc::er_no_such_table);
    BOOST_TEST(res[1].diag() == create_server_diag("Table 'boost_mysql_integtests.bad_table' doesn't exist"));
    BOOST_TEST(res[2].error() == common_server_errc::er_empty_query);
    BOOST_TEST(res[2].diag() == create_server_diag("Query was empty"));
    BOOST_TEST(res[3].as_results().rows() == makerows(1, 42), per_element());
    BOOST_TEST(res[4].error() == common_server_errc::er_unknown_character_set);
    BOOST_TEST(res[4].diag() == create_server_diag("Unknown character set: 'bad_charset'"));
    BOOST_TEST(res[5].as_results().rows() == makerows(1, "abc"), per_element());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace