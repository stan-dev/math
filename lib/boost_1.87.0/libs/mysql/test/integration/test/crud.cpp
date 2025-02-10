//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/results.hpp>

#include <boost/test/unit_test_suite.hpp>

#include "test_common/create_basic.hpp"
#include "test_common/network_result.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/metadata_validator.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::test_tools::per_element;

namespace {

BOOST_AUTO_TEST_SUITE(test_crud)

// Other SELECT statements are already covered
BOOST_FIXTURE_TEST_CASE(query_empty_select, any_connection_fixture)
{
    connect();

    // Issue query
    results result;
    conn.async_execute("SELECT * FROM empty_table", result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST(result.size() == 1u);
    validate_2fields_meta(result.meta(), "empty_table");
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(query_empty_select_multifn, any_connection_fixture)
{
    connect();

    // Issue query
    execution_state st;
    conn.async_start_execution("SELECT * FROM empty_table", st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    validate_2fields_meta(st.meta(), "empty_table");

    // Read eof
    auto rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rv == rows_view(), per_element());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
}

BOOST_FIXTURE_TEST_CASE(query_insert, any_connection_fixture)
{
    connect();
    start_transaction();

    // Issue query
    constexpr const char*
        query = "INSERT INTO inserts_table (field_varchar, field_date) VALUES ('v0', '2010-10-11')";
    results result;
    conn.async_execute(query, result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST(result.size() == 1u);
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.affected_rows() == 1u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() > 0u);
    BOOST_TEST(result.info() == "");

    // Verify insertion took place
    conn.async_execute("SELECT COUNT(*) FROM inserts_table", result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, 1), per_element());
}

BOOST_FIXTURE_TEST_CASE(query_update, any_connection_fixture)
{
    connect();
    start_transaction();

    // Issue the query
    results result;
    conn.async_execute("UPDATE updates_table SET field_int = field_int+10", result, as_netresult)
        .validate_no_error();

    // Validate results
    BOOST_TEST(result.size() == 1u);
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.affected_rows() == 2u);  // there are 3 rows, but 1 has field_int = NULL
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "Rows matched: 3  Changed: 2  Warnings: 0");

    // Validate it took effect
    conn.async_execute("SELECT field_int FROM updates_table WHERE field_varchar = 'f0'", result, as_netresult)
        .validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, 52), per_element());  // initial value was 42
}

BOOST_FIXTURE_TEST_CASE(query_delete, any_connection_fixture)
{
    connect();
    start_transaction();

    // Issue the query
    results result;
    conn.async_execute("DELETE FROM updates_table", result, as_netresult).validate_no_error();

    // Validate results
    BOOST_TEST(result.size() == 1u);
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.affected_rows() == 3u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Validate it took effect
    conn.async_execute("SELECT COUNT(*) FROM updates_table", result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, 0), per_element());
}

BOOST_FIXTURE_TEST_CASE(statement_update, any_connection_fixture)
{
    connect();
    start_transaction();

    // Prepare the statement
    constexpr const char* sql = "UPDATE updates_table SET field_int = ? WHERE field_varchar = ?";
    auto stmt = conn.async_prepare_statement(sql, as_netresult).get();
    BOOST_TEST(stmt.num_params() == 2u);

    // Execute it
    results result;
    conn.execute(stmt.bind(200, "f0"), result);
    BOOST_TEST(result.size() == 1u);
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.affected_rows() == 1u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "Rows matched: 1  Changed: 1  Warnings: 0");

    // Verify that it took effect
    conn.async_execute("SELECT field_int FROM updates_table WHERE field_varchar = 'f0'", result, as_netresult)
        .validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, 200), per_element());

    // Close the statement
    conn.async_close_statement(stmt, as_netresult).validate_no_error();
}

BOOST_FIXTURE_TEST_CASE(statement_delete, any_connection_fixture)
{
    connect();
    start_transaction();

    // Prepare the statement
    constexpr const char* sql = "DELETE FROM updates_table WHERE field_varchar = ?";
    auto stmt = conn.async_prepare_statement(sql, as_netresult).get();
    BOOST_TEST(stmt.num_params() == 1u);

    // Execute it
    results result;
    conn.async_execute(stmt.bind("f0"), result, as_netresult).validate_no_error();
    BOOST_TEST(result.size() == 1u);
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.affected_rows() == 1u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Validate it took effect
    conn.async_execute("SELECT COUNT(*) FROM updates_table", result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, 2), per_element());
}

BOOST_AUTO_TEST_SUITE_END()  // test_crud

}  // namespace
