//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/create_basic.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/metadata_validator.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::test_tools::per_element;

namespace {

BOOST_AUTO_TEST_SUITE(test_prepared_statements)

BOOST_FIXTURE_TEST_CASE(multiple_executions, any_connection_fixture)
{
    connect();
    results result;

    // Prepare a statement
    constexpr const char* sql = "SELECT * FROM two_rows_table WHERE id = ? OR field_varchar = ?";
    auto stmt = conn.async_prepare_statement(sql, as_netresult).get();

    // Execute it. Only one row will be returned (because of the id)
    conn.async_execute(stmt.bind(1, "non_existent"), result, as_netresult).validate_no_error();
    validate_2fields_meta(result.meta(), "two_rows_table");
    BOOST_TEST(result.rows() == makerows(2, 1, "f0"), per_element());

    // Execute it again, but with different values. This time, two rows are returned
    conn.async_execute(stmt.bind(1, "f1"), result, as_netresult).validate_no_error();
    validate_2fields_meta(result.meta(), "two_rows_table");
    BOOST_TEST(result.rows() == makerows(2, 1, "f0", 2, "f1"), per_element());

    // Close it
    conn.async_close_statement(stmt, as_netresult).validate_no_error();
}

BOOST_FIXTURE_TEST_CASE(multiple_statements, any_connection_fixture)
{
    connect();
    start_transaction();
    results result;

    // Prepare an update and a select
    constexpr const char* update_sql = "UPDATE updates_table SET field_int = ? WHERE field_varchar = ?";
    constexpr const char* select_sql = "SELECT field_int FROM updates_table WHERE field_varchar = ?";
    auto stmt_update = conn.async_prepare_statement(update_sql, as_netresult).get();
    auto stmt_select = conn.async_prepare_statement(select_sql, as_netresult).get();
    BOOST_TEST(stmt_update.num_params() == 2u);
    BOOST_TEST(stmt_select.num_params() == 1u);
    BOOST_TEST(stmt_update.id() != stmt_select.id());

    // Execute update
    conn.async_execute(stmt_update.bind(210, "f0"), result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 1u);

    // Execute select
    conn.async_execute(stmt_select.bind("f0"), result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, 210), per_element());

    // Execute update again
    conn.async_execute(stmt_update.bind(220, "f0"), result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 1u);

    // Update no longer needed, close it
    conn.async_close_statement(stmt_update, as_netresult).validate_no_error();

    // Execute select again
    conn.async_execute(stmt_select.bind("f0"), result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, 220), per_element());

    // Close select
    conn.async_close_statement(stmt_select, as_netresult).validate_no_error();
}

BOOST_FIXTURE_TEST_CASE(statement_without_params, any_connection_fixture)
{
    connect();

    // Prepare the statement
    auto stmt = conn.async_prepare_statement("SELECT * FROM empty_table", as_netresult).get();
    BOOST_TEST(stmt.valid());
    BOOST_TEST(stmt.num_params() == 0u);

    // Execute doesn't error
    results result;
    conn.async_execute(stmt.bind(), result, as_netresult).validate_no_error();
    validate_2fields_meta(result.meta(), "empty_table");
    BOOST_TEST(result.rows() == rows(), per_element());
}

BOOST_FIXTURE_TEST_CASE(wrong_num_params, any_connection_fixture)
{
    connect();

    // Prepare the statement
    auto stmt = conn.async_prepare_statement("SELECT * FROM empty_table", as_netresult).get();
    BOOST_TEST(stmt.valid());
    BOOST_TEST(stmt.num_params() == 0u);

    // Execute fails appropriately
    results result;
    conn.async_execute(stmt.bind(42), result, as_netresult).validate_error(client_errc::wrong_num_params);
}

// Note: multifn query is already covered in spotchecks
BOOST_FIXTURE_TEST_CASE(multifn, any_connection_fixture)
{
    connect();

    // Prepare the statement
    auto stmt = conn.async_prepare_statement("SELECT * FROM three_rows_table", as_netresult).get();

    // Execute it
    execution_state st;
    conn.async_start_execution(stmt.bind(), st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // We don't know how many rows there will be in each batch,
    // but they will come in order
    std::size_t call_count = 0;
    std::vector<row> all_rows;
    while (st.should_read_rows() && call_count <= 4)
    {
        ++call_count;
        for (row_view rv : conn.async_read_some_rows(st, as_netresult).get())
            all_rows.emplace_back(rv);
    }

    // Verify rows
    BOOST_TEST(all_rows == makerows(2, 1, "f0", 2, "f1", 3, "f2"), per_element());

    // Verify eof
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");

    // Reading again does nothing
    auto rws = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rws == rows(), per_element());
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
}

BOOST_AUTO_TEST_SUITE_END()  // test_prepared_statements

}  // namespace
