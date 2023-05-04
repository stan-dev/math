//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/tcp.hpp>

#include <boost/test/unit_test_suite.hpp>

#include "integration_test_common.hpp"
#include "tcp_network_fixture.hpp"
#include "test_common.hpp"

using namespace boost::mysql::test;
using boost::mysql::results;

namespace {

BOOST_AUTO_TEST_SUITE(test_crud)

// Other SELECT statements are already covered
BOOST_FIXTURE_TEST_CASE(query_empty_select, tcp_network_fixture)
{
    connect();

    // Issue query
    results result;
    conn.query("SELECT * FROM empty_table", result);

    // Verify results
    validate_2fields_meta(result.meta(), "empty_table");
    BOOST_TEST(result.rows().empty());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(query_insert, tcp_network_fixture)
{
    connect();
    start_transaction();

    // Issue query
    results result;
    conn.query("INSERT INTO inserts_table (field_varchar, field_date) VALUES ('v0', '2010-10-11')", result);

    // Verify results
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows().empty());
    BOOST_TEST(result.affected_rows() == 1u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() > 0u);
    BOOST_TEST(result.info() == "");

    // Verify insertion took place
    conn.query("SELECT COUNT(*) FROM inserts_table", result);
    BOOST_TEST(result.rows().at(0).at(0).as_int64() == 1);
}

BOOST_FIXTURE_TEST_CASE(query_update, tcp_network_fixture)
{
    connect();
    start_transaction();

    // Issue the query
    results result;
    conn.query("UPDATE updates_table SET field_int = field_int+10", result);

    // Validate results
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows().empty());
    BOOST_TEST(result.affected_rows() == 2u);  // there are 3 rows, but 1 has field_int = NULL
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "Rows matched: 3  Changed: 2  Warnings: 0");

    // Validate it took effect
    conn.query("SELECT field_int FROM updates_table WHERE field_varchar = 'f0'", result);
    BOOST_TEST(result.rows().at(0).at(0).as_int64() == 52);  // initial value was 42
}

BOOST_FIXTURE_TEST_CASE(query_delete, tcp_network_fixture)
{
    connect();
    start_transaction();

    // Issue the query
    results result;
    conn.query("DELETE FROM updates_table", result);

    // Validate results
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows().empty());
    BOOST_TEST(result.affected_rows() == 3u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Validate it took effect
    conn.query("SELECT COUNT(*) FROM updates_table", result);
    BOOST_TEST(result.rows().at(0).at(0).as_int64() == 0);
}

BOOST_FIXTURE_TEST_CASE(statement_update, tcp_network_fixture)
{
    connect();
    start_transaction();

    // Prepare the statement
    auto stmt = conn.prepare_statement("UPDATE updates_table SET field_int = ? WHERE field_varchar = ?");
    BOOST_TEST(stmt.num_params() == 2u);

    // Execute it
    results result;
    conn.execute_statement(stmt, std::make_tuple(200, "f0"), result);
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows().empty());
    BOOST_TEST(result.affected_rows() == 1u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "Rows matched: 1  Changed: 1  Warnings: 0");

    // Verify that it took effect
    conn.query("SELECT field_int FROM updates_table WHERE field_varchar = 'f0'", result);
    BOOST_TEST(result.rows().at(0).at(0).as_int64() == 200);

    // Close the statement
    conn.close_statement(stmt);
}

BOOST_FIXTURE_TEST_CASE(statement_delete, tcp_network_fixture)
{
    connect();
    start_transaction();

    // Prepare the statement
    auto stmt = conn.prepare_statement("DELETE FROM updates_table WHERE field_varchar = ?");
    BOOST_TEST(stmt.num_params() == 1u);

    // Execute it
    results result;
    conn.execute_statement(stmt, std::make_tuple("f0"), result);
    BOOST_TEST(result.meta().empty());
    BOOST_TEST(result.rows().empty());
    BOOST_TEST(result.affected_rows() == 1u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Validate it took effect
    conn.query("SELECT COUNT(*) FROM updates_table", result);
    BOOST_TEST(result.rows().at(0).at(0).as_int64() == 2);
}

BOOST_AUTO_TEST_SUITE_END()  // test_crud

}  // namespace
