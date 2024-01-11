//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/tcp.hpp>

#include <boost/test/unit_test.hpp>

#include <tuple>

#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_integration/common.hpp"
#include "test_integration/er_connection.hpp"
#include "test_integration/tcp_network_fixture.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_SUITE(test_prepared_statements)

BOOST_FIXTURE_TEST_CASE(multiple_executions, tcp_network_fixture)
{
    connect();

    // Prepare a statement
    auto stmt = conn.prepare_statement("SELECT * FROM two_rows_table WHERE id = ? OR field_varchar = ?");

    // Execute it. Only one row will be returned (because of the id)
    results result;
    conn.execute(stmt.bind(1, "non_existent"), result);
    validate_2fields_meta(result.meta(), "two_rows_table");
    BOOST_TEST_REQUIRE(result.rows().size() == 1u);
    BOOST_TEST((result.rows()[0] == makerow(1, "f0")));

    // Execute it again, but with different values. This time, two rows are returned
    conn.execute(stmt.bind(1, "f1"), result);
    validate_2fields_meta(result.meta(), "two_rows_table");
    BOOST_TEST_REQUIRE(result.rows().size() == 2u);
    BOOST_TEST((result.rows()[0] == makerow(1, "f0")));
    BOOST_TEST((result.rows()[1] == makerow(2, "f1")));

    // Close it
    conn.close_statement(stmt);
}

BOOST_FIXTURE_TEST_CASE(multiple_statements, tcp_network_fixture)
{
    connect();
    start_transaction();

    // Prepare an update and a select
    results result;
    auto stmt_update = conn.prepare_statement("UPDATE updates_table SET field_int = ? WHERE field_varchar = ?"
    );
    auto stmt_select = conn.prepare_statement("SELECT field_int FROM updates_table WHERE field_varchar = ?");
    BOOST_TEST(stmt_update.num_params() == 2u);
    BOOST_TEST(stmt_select.num_params() == 1u);
    BOOST_TEST(stmt_update.id() != stmt_select.id());

    // Execute update
    conn.execute(stmt_update.bind(210, "f0"), result);
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 1u);

    // Execute select
    conn.execute(stmt_select.bind("f0"), result);
    BOOST_TEST(result.rows().size() == 1u);
    BOOST_TEST(result.rows()[0] == makerow(210));

    // Execute update again
    conn.execute(stmt_update.bind(220, "f0"), result);
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 1u);

    // Update no longer needed, close it
    conn.close_statement(stmt_update);

    // Execute select again
    conn.execute(stmt_select.bind("f0"), result);
    BOOST_TEST(result.rows().size() == 1u);
    BOOST_TEST(result.rows()[0] == makerow(220));

    // Close select
    conn.close_statement(stmt_select);
}

BOOST_FIXTURE_TEST_CASE(statement_without_params, tcp_network_fixture)
{
    connect();

    // Prepare the statement
    auto stmt = conn.prepare_statement("SELECT * FROM empty_table");
    BOOST_TEST(stmt.valid());
    BOOST_TEST(stmt.num_params() == 0u);

    // Execute doesn't error
    results result;
    conn.execute(stmt.bind(), result);
    validate_2fields_meta(result.meta(), "empty_table");
    BOOST_TEST(result.rows().size() == 0u);
}

BOOST_FIXTURE_TEST_CASE(wrong_num_params, tcp_network_fixture)
{
    connect();

    // Prepare the statement
    auto stmt = conn.prepare_statement("SELECT * FROM empty_table");
    BOOST_TEST(stmt.valid());
    BOOST_TEST(stmt.num_params() == 0u);

    // Execute fails appropriately
    error_code ec;
    diagnostics diag;
    results result;
    conn.execute(stmt.bind(42), result, ec, diag);
    BOOST_TEST(ec == client_errc::wrong_num_params);
}

// Note: multifn query is already covered in spotchecks
BOOST_FIXTURE_TEST_CASE(multifn, tcp_network_fixture)
{
    connect();

    // Prepare the statement
    auto stmt = conn.prepare_statement("SELECT * FROM three_rows_table");

    // Execute it
    execution_state st;
    conn.start_execution(stmt.bind(), st);
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // We don't know how many rows there will be in each batch,
    // but they will come in order
    std::size_t call_count = 0;
    std::vector<row> all_rows;
    while (st.should_read_rows() && call_count <= 4)
    {
        ++call_count;
        for (row_view rv : conn.read_some_rows(st))
            all_rows.emplace_back(rv);
    }

    // Verify rows
    BOOST_TEST_REQUIRE(all_rows.size() == 3u);
    BOOST_TEST(all_rows[0] == makerow(1, "f0"));
    BOOST_TEST(all_rows[1] == makerow(2, "f1"));
    BOOST_TEST(all_rows[2] == makerow(3, "f2"));

    // Verify eof
    validate_eof(st);

    // Reading again does nothing
    auto rows = conn.read_some_rows(st);
    BOOST_TEST(rows.empty());
    validate_eof(st);
}

BOOST_AUTO_TEST_SUITE_END()  // test_prepared_statements

}  // namespace
