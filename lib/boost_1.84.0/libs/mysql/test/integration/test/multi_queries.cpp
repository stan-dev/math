//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/results.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_common/validate_string_contains.hpp"
#include "test_integration/common.hpp"
#include "test_integration/tcp_network_fixture.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_SUITE(test_multi_queries)

BOOST_FIXTURE_TEST_CASE(empty_results, tcp_network_fixture)
{
    params.set_multi_queries(true);
    connect();
    start_transaction();

    results result;
    conn.execute(
        R"%(
            INSERT INTO inserts_table (field_varchar) VALUES ('abc');
            INSERT INTO inserts_table (field_varchar) VALUES ('def');
            DELETE FROM updates_table
        )%",
        result
    );

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

BOOST_FIXTURE_TEST_CASE(data_results, tcp_network_fixture)
{
    params.set_multi_queries(true);
    connect();
    start_transaction();

    results result;
    conn.execute(
        R"%(
            SELECT * FROM one_row_table;
            SELECT * FROM empty_table;
            DELETE FROM updates_table
        )%",
        result
    );

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

BOOST_FIXTURE_TEST_CASE(error_not_enabled, tcp_network_fixture)
{
    connect();

    results result;
    error_code err;
    diagnostics diag;

    conn.execute("SELECT 1; SELECT 2", result, err, diag);
    BOOST_TEST(err == common_server_errc::er_parse_error);
    validate_string_contains(diag.server_message(), {"you have an error in your sql syntax"});
}

BOOST_AUTO_TEST_SUITE_END()  // test_crud

}  // namespace
