//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/results.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/create_basic.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/metadata_validator.hpp"
#include "test_integration/tcp_connection_fixture.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::test_tools::per_element;

namespace {

BOOST_AUTO_TEST_SUITE(test_multi_queries)

BOOST_FIXTURE_TEST_CASE(empty_results, any_connection_fixture)
{
    // Setup
    constexpr const char* query =
        "INSERT INTO inserts_table (field_varchar) VALUES ('abc');"
        "INSERT INTO inserts_table (field_varchar) VALUES ('def');"
        "DELETE FROM updates_table";
    connect(connect_params_builder().disable_ssl().multi_queries(true).build());
    start_transaction();

    // Run the query
    results result;
    conn.async_execute(query, result, as_netresult).validate_no_error();

    // Validate results
    BOOST_TEST_REQUIRE(result.size() == 3u);

    BOOST_TEST(result[0].meta().size() == 0u);
    BOOST_TEST(result[0].rows() == rows());
    BOOST_TEST(result[0].affected_rows() == 1u);
    BOOST_TEST(result[0].warning_count() == 0u);
    BOOST_TEST(result[0].last_insert_id() > 0u);
    BOOST_TEST(result[0].info() == "");
    BOOST_TEST(!result[0].is_out_params());

    BOOST_TEST(result[1].meta().size() == 0u);
    BOOST_TEST(result[1].rows() == rows());
    BOOST_TEST(result[1].affected_rows() == 1u);
    BOOST_TEST(result[1].warning_count() == 0u);
    BOOST_TEST(result[1].last_insert_id() > 0u);
    BOOST_TEST(result[1].info() == "");
    BOOST_TEST(!result[1].is_out_params());

    BOOST_TEST(result[2].meta().size() == 0u);
    BOOST_TEST(result[2].rows() == rows());
    BOOST_TEST(result[2].affected_rows() == 3u);
    BOOST_TEST(result[2].warning_count() == 0u);
    BOOST_TEST(result[2].last_insert_id() == 0u);
    BOOST_TEST(result[2].info() == "");
    BOOST_TEST(!result[2].is_out_params());
}

BOOST_FIXTURE_TEST_CASE(data_results, any_connection_fixture)
{
    // Setup
    connect(connect_params_builder().disable_ssl().multi_queries(true).build());
    start_transaction();
    constexpr const char* query =
        "SELECT * FROM one_row_table;"
        "SELECT * FROM empty_table;"
        "DELETE FROM updates_table";

    // Execute
    results result;
    conn.async_execute(query, result, as_netresult).validate_no_error();

    // Validate results
    BOOST_TEST_REQUIRE(result.size() == 3u);

    validate_2fields_meta(result[0].meta(), "one_row_table");
    BOOST_TEST(result[0].rows() == makerows(2, 1, "f0"));
    BOOST_TEST(result[0].affected_rows() == 0u);
    BOOST_TEST(result[0].warning_count() == 0u);
    BOOST_TEST(result[0].last_insert_id() == 0u);
    BOOST_TEST(result[0].info() == "");
    BOOST_TEST(!result[0].is_out_params());

    validate_2fields_meta(result[1].meta(), "empty_table");
    BOOST_TEST(result[1].rows() == makerows(2));
    BOOST_TEST(result[1].affected_rows() == 0u);
    BOOST_TEST(result[1].warning_count() == 0u);
    BOOST_TEST(result[1].last_insert_id() == 0u);
    BOOST_TEST(result[1].info() == "");
    BOOST_TEST(!result[1].is_out_params());

    BOOST_TEST(result[2].meta().size() == 0u);
    BOOST_TEST(result[2].rows() == rows());
    BOOST_TEST(result[2].affected_rows() == 3u);
    BOOST_TEST(result[2].warning_count() == 0u);
    BOOST_TEST(result[2].last_insert_id() == 0u);
    BOOST_TEST(result[2].info() == "");
    BOOST_TEST(!result[2].is_out_params());
}

BOOST_FIXTURE_TEST_CASE(error_not_enabled, any_connection_fixture)
{
    // Setup
    connect(connect_params_builder().disable_ssl().build());

    // Execute fails
    results result;
    conn.async_execute("SELECT 1; SELECT 2", result, as_netresult)
        .validate_error_contains(
            common_server_errc::er_parse_error,
            {"you have an error in your sql syntax"}
        );
}

// Spotcheck: old connection can do multi-queries
BOOST_FIXTURE_TEST_CASE(tcp_connection_enable, tcp_connection_fixture)
{
    connect(connect_params_builder().multi_queries(true).build_hparams());

    // Execute succeeds
    results result;
    conn.async_execute("SELECT 1; SELECT 2", result, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(result.size() == 2u);
    BOOST_TEST(result.at(0).rows() == makerows(1, 1), per_element());
    BOOST_TEST(result.at(1).rows() == makerows(1, 2), per_element());
}

BOOST_FIXTURE_TEST_CASE(tcp_connection_disabled, tcp_connection_fixture)
{
    // Setup
    connect();

    // Execute succeeds
    results result;
    conn.async_execute("SELECT 1; SELECT 2", result, as_netresult)
        .validate_error_contains(
            common_server_errc::er_parse_error,
            {"you have an error in your sql syntax"}
        );
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
