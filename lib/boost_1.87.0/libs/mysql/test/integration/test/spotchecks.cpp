//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/pipeline.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/tcp.hpp>

#include <boost/core/span.hpp>
#include <boost/optional/optional.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>

#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/metadata_validator.hpp"
#include "test_integration/spotchecks_helpers.hpp"
#include "test_integration/static_rows.hpp"
#include "test_integration/tcp_connection_fixture.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
namespace asio = boost::asio;
using boost::test_tools::per_element;

// These tests aim to cover the 4 overloads we have for each network function,
// both for the old templated connection and for any_connection.
// A success and an error case is included for each function.

namespace {

BOOST_AUTO_TEST_SUITE(test_spotchecks)

// Functions for each connection type
static auto fns_connection = network_functions_connection::all();
static auto fns_any = network_functions_any::all();

// Simplify defining test cases involving both connection types
// The resulting test case gets the fixture (fix) and the network functions (fn) as args
#define BOOST_MYSQL_SPOTCHECK_TEST(name)                    \
    template <class FixtureType>                            \
    void do_##name(FixtureType& fix);                       \
    BOOST_DATA_TEST_CASE(name##_connection, fns_connection) \
    {                                                       \
        netfn_fixture_connection fix{sample};               \
        do_##name(fix);                                     \
    }                                                       \
    BOOST_DATA_TEST_CASE(name##_any, fns_any)               \
    {                                                       \
        netfn_fixture_any fix{sample};                      \
        do_##name(fix);                                     \
    }                                                       \
    template <class FixtureType>                            \
    void do_##name(FixtureType& fix)

// prepare statement
BOOST_MYSQL_SPOTCHECK_TEST(prepare_statement_success)
{
    // Setup
    fix.connect();

    // Call the function
    statement stmt = fix.net.prepare_statement(fix.conn, "SELECT * FROM empty_table WHERE id IN (?, ?)")
                         .get();

    // Validate the result
    BOOST_TEST_REQUIRE(stmt.valid());
    BOOST_TEST(stmt.id() > 0u);
    BOOST_TEST(stmt.num_params() == 2u);

    // It can be executed
    results result;
    fix.net.execute_statement(fix.conn, stmt.bind(10, 20), result).validate_no_error();
    BOOST_TEST(result.rows().empty());
}

BOOST_MYSQL_SPOTCHECK_TEST(prepare_statement_error)
{
    // Setup
    fix.connect();

    // Call the function
    fix.net.prepare_statement(fix.conn, "SELECT * FROM bad_table WHERE id IN (?, ?)")
        .validate_error(
            common_server_errc::er_no_such_table,
            "Table 'boost_mysql_integtests.bad_table' doesn't exist"
        );
}

// execute
BOOST_MYSQL_SPOTCHECK_TEST(execute_success)
{
    // Setup
    fix.connect();

    // Call the function
    results result;
    fix.net.execute_query(fix.conn, "SELECT 'hello', 42", result).validate_no_error();

    // Check results
    BOOST_TEST(result.rows().size() == 1u);
    BOOST_TEST(result.rows()[0] == makerow("hello", 42), per_element());
    BOOST_TEST(result.meta().size() == 2u);
}

BOOST_MYSQL_SPOTCHECK_TEST(execute_error)
{
    // Setup
    fix.connect();

    // Call the function
    results result;
    fix.net.execute_query(fix.conn, "SELECT field_varchar, field_bad FROM one_row_table", result)
        .validate_error_contains(common_server_errc::er_bad_field_error, {"unknown column", "field_bad"});
}

// Start execution
BOOST_MYSQL_SPOTCHECK_TEST(start_execution_success)
{
    // Setup
    fix.connect();

    // Call the function
    execution_state st;
    fix.net.start_execution(fix.conn, "SELECT * FROM empty_table", st).validate_no_error();

    // Check results
    BOOST_TEST(st.should_read_rows());
    validate_2fields_meta(st.meta(), "empty_table");
}

BOOST_MYSQL_SPOTCHECK_TEST(start_execution_error)
{
    // Setup
    fix.connect();

    // Call the function
    execution_state st;
    fix.net.start_execution(fix.conn, "SELECT field_varchar, field_bad FROM one_row_table", st)
        .validate_error_contains(common_server_errc::er_bad_field_error, {"unknown column", "field_bad"});
}

// Close statement
BOOST_MYSQL_SPOTCHECK_TEST(close_statement_success)
{
    // Setup
    fix.connect();

    // Prepare a statement
    statement stmt = fix.net.prepare_statement(fix.conn, "SELECT * FROM empty_table WHERE id IN (?, ?)")
                         .get();

    // Close the statement
    fix.net.close_statement(fix.conn, stmt).validate_no_error();

    // The statement is no longer valid
    results result;
    fix.net.execute_statement(fix.conn, stmt.bind(1, 2), result).validate_any_error();
}

BOOST_MYSQL_SPOTCHECK_TEST(close_statement_error)
{
    // Setup
    fix.connect();

    // Prepare a statement
    statement stmt = fix.net.prepare_statement(fix.conn, "SELECT * FROM empty_table WHERE id IN (?, ?)")
                         .get();

    // Close the connection
    fix.net.close(fix.conn).validate_no_error();

    // Close the statement fails, as it requires communication with the server
    fix.net.close_statement(fix.conn, stmt).validate_any_error();
}

// Read resultset head
BOOST_MYSQL_SPOTCHECK_TEST(read_resultset_head_success)
{
    // Setup
    fix.connect(connect_params_builder().multi_queries(true));

    // Generate an execution state
    execution_state st;
    fix.net.start_execution(fix.conn, "SELECT 4.2e0; SELECT * FROM empty_table", st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // Read the 1st resultset
    auto rws = fix.net.read_some_rows(fix.conn, st).get();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(st.meta()[0].type() == column_type::double_);
    BOOST_TEST(rws == makerows(1, 4.2), per_element());

    // Read head
    fix.net.read_resultset_head(fix.conn, st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    validate_2fields_meta(st.meta(), "empty_table");

    // Reading head again does nothing
    fix.net.read_resultset_head(fix.conn, st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    validate_2fields_meta(st.meta(), "empty_table");

    // We can read rows now
    rws = fix.net.read_some_rows(fix.conn, st).get();
    BOOST_TEST(rws == rows(), per_element());
}

BOOST_MYSQL_SPOTCHECK_TEST(read_resultset_head_error)
{
    // Setup
    fix.connect(connect_params_builder().multi_queries(true));

    // Generate an execution state
    execution_state st;
    fix.net.start_execution(fix.conn, "SELECT * FROM empty_table; SELECT bad_field FROM one_row_table", st)
        .validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // Read the OK packet to finish 1st resultset
    fix.net.read_some_rows(fix.conn, st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_head());

    // Read head for the 2nd resultset. This one contains an error, which is detected when reading head.
    fix.net.read_resultset_head(fix.conn, st)
        .validate_error(common_server_errc::er_bad_field_error, "Unknown column 'bad_field' in 'field list'");
}

// Read some rows. No error spotcheck here
// TODO: write a unit test with this
BOOST_MYSQL_SPOTCHECK_TEST(read_some_rows_success)
{
    // Setup
    fix.connect();

    // Generate an execution state
    execution_state st;
    fix.net.start_execution(fix.conn, "SELECT * FROM one_row_table", st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // Read once. st may or may not be complete, depending
    // on how the buffer reallocated memory
    auto rows = fix.net.read_some_rows(fix.conn, st).get();
    BOOST_TEST((rows == makerows(2, 1, "f0")));

    // Reading again should complete st
    rows = fix.net.read_some_rows(fix.conn, st).get();
    BOOST_TEST(rows.empty());
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");

    // Reading again does nothing
    rows = fix.net.read_some_rows(fix.conn, st).get();
    BOOST_TEST(rows.empty());
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
}

// Ping
BOOST_MYSQL_SPOTCHECK_TEST(ping_success)
{
    // Setup
    fix.connect();

    // Success
    fix.net.ping(fix.conn).validate_no_error();
}

BOOST_MYSQL_SPOTCHECK_TEST(ping_error)
{
    // Pinging an unconnected connection fails
    fix.net.ping(fix.conn).validate_any_error();
}

// Reset connection
BOOST_MYSQL_SPOTCHECK_TEST(reset_connection_success)
{
    // Setup
    fix.connect();

    // Set some variable
    results result;
    fix.net.execute_query(fix.conn, "SET @myvar = 42", result).validate_no_error();

    // Reset connection
    fix.net.reset_connection(fix.conn).validate_no_error();

    // The variable has been reset
    fix.net.execute_query(fix.conn, "SELECT @myvar", result).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, nullptr), per_element());
}

BOOST_MYSQL_SPOTCHECK_TEST(reset_connection_error)
{
    // Resetting an unconnected connection fails
    fix.net.reset_connection(fix.conn).validate_any_error();
}

// Close connection
// TODO: make a unit test with double close and an error case
BOOST_MYSQL_SPOTCHECK_TEST(close_success)
{
    // Setup
    fix.connect();

    // Close
    fix.net.close(fix.conn).validate_no_error();

    // We are no longer able to query
    results result;
    fix.net.execute_query(fix.conn, "SELECT 1", result).validate_any_error();

    // Closing again returns OK (and does nothing)
    fix.net.close(fix.conn).validate_no_error();
}

#ifdef BOOST_MYSQL_CXX14
// Execute (static) - errors are already covered by the other tests
BOOST_MYSQL_SPOTCHECK_TEST(execute_static_success)
{
    // Setup
    fix.connect();
    static_results_t result;

    // Execute the function
    fix.net.execute_static(fix.conn, "CALL sp_spotchecks()", result).validate_no_error();

    // Check
    row_multifield expected[]{
        {boost::optional<float>(1.1f), 11, std::string("aaa")}
    };
    BOOST_TEST(result.rows<0>() == expected, per_element());
}

// start_execution, read_resultset_head, read_some_rows success
BOOST_MYSQL_SPOTCHECK_TEST(start_execution_and_followups_static_success)
{
    // Setup
    fix.connect();
    static_state_t st;

    // Start
    fix.net.start_execution_static(fix.conn, "CALL sp_spotchecks()", st).validate_no_error();
    BOOST_TEST(st.should_read_rows());

    // Read r1 rows
    std::array<row_multifield, 2> storage;
    std::size_t num_rows = fix.net.read_some_rows_static_1(fix.conn, st, storage).get();
    row_multifield expected_multifield{boost::optional<float>(1.1f), 11, std::string("aaa")};  // MSVC 14.1
    BOOST_TEST(num_rows == 1u);
    BOOST_TEST(storage[0] == expected_multifield);

    // Ensure we're in the next resultset
    num_rows = fix.net.read_some_rows_static_1(fix.conn, st, storage).get();
    BOOST_TEST(num_rows == 0u);
    BOOST_TEST(st.should_read_head());

    // Read r2 head
    fix.net.read_resultset_head_static(fix.conn, st).validate_no_error();
    BOOST_TEST(st.should_read_rows());

    // Read r2 rows
    std::array<row_2fields, 2> storage2;
    num_rows = fix.net.read_some_rows_static_2(fix.conn, st, storage2).get();
    BOOST_TEST(num_rows == 1u);
    row_2fields expected_2fields{1, std::string("f0")};
    BOOST_TEST(storage2[0] == expected_2fields);

    // Ensure we're in the next resultset
    num_rows = fix.net.read_some_rows_static_2(fix.conn, st, storage2).get();
    BOOST_TEST(num_rows == 0u);
    BOOST_TEST(st.should_read_head());

    // Read r3 head (empty)
    fix.net.read_resultset_head_static(fix.conn, st).validate_no_error();
    BOOST_TEST(st.complete());
}

// read_some_rows failure (the other error cases are already widely tested)
BOOST_MYSQL_SPOTCHECK_TEST(read_some_rows_static_error)
{
    // Setup
    fix.connect();
    static_state_t st;

    // Start
    fix.net.start_execution_static(fix.conn, "SELECT * FROM multifield_table WHERE id = 42", st)
        .validate_no_error();
    BOOST_TEST(st.should_read_rows());

    // No rows matched, so reading rows reads the OK packet. This will report the num resultsets mismatch
    std::array<row_multifield, 2> storage;
    fix.net.read_some_rows_static_1(fix.conn, st, storage)
        .validate_error(client_errc::num_resultsets_mismatch);
}
#endif

//
// functions specific to old connection
//

// Handshake
BOOST_DATA_TEST_CASE_F(tcp_connection_fixture, handshake_success, fns_connection)
{
    // Setup
    const network_functions_connection& fn = sample;
    fn.connect_stream(conn.stream(), get_tcp_endpoint()).validate_no_error_nodiag();

    // Call the function
    fn.handshake(conn, connect_params_builder().build_hparams()).validate_no_error();

    // The connection is usable
    fn.ping(conn).validate_no_error();
}

BOOST_DATA_TEST_CASE_F(tcp_connection_fixture, handshake_error, fns_connection)
{
    // Setup
    const network_functions_connection& fn = sample;
    fn.connect_stream(conn.stream(), get_tcp_endpoint()).validate_no_error_nodiag();

    // Call the function
    fn.handshake(conn, connect_params_builder().database("bad_db").build_hparams())
        .validate_error(
            common_server_errc::er_dbaccess_denied_error,
            "Access denied for user 'integ_user'@'%' to database 'bad_db'"
        );
}

// Connect
BOOST_DATA_TEST_CASE_F(tcp_connection_fixture, connect_connection_success, fns_connection)
{
    // Setup
    const network_functions_connection& fn = sample;

    // Call the function
    fn.connect(conn, get_tcp_endpoint(), connect_params_builder().build_hparams()).validate_no_error();

    // The connection is usable
    fn.ping(conn).validate_no_error();
}

BOOST_DATA_TEST_CASE_F(tcp_connection_fixture, connect_connection_error, fns_connection)
{
    // Setup
    const network_functions_connection& fn = sample;

    // Call the function
    fn.connect(conn, get_tcp_endpoint(), connect_params_builder().database("bad_db").build_hparams())
        .validate_error(
            common_server_errc::er_dbaccess_denied_error,
            "Access denied for user 'integ_user'@'%' to database 'bad_db'"
        );
}

// Quit. TODO: write a unit test spotcheck for errors
BOOST_DATA_TEST_CASE_F(tcp_connection_fixture, quit_success, fns_connection)
{
    // Setup
    const network_functions_connection& fn = sample;
    connect();

    // Quit
    fn.quit(conn).validate_no_error();
}

//
// functions specific to any_connection
//

// Connect
BOOST_DATA_TEST_CASE_F(any_connection_fixture, connect_any_success, fns_any)
{
    // Setup
    const network_functions_any& fn = sample;

    // Call the function
    fn.connect(conn, connect_params_builder().ssl(ssl_mode::require).build()).validate_no_error();

    // The connection is usable
    fn.ping(conn).validate_no_error();

    // Closing works
    fn.close(conn).validate_no_error();
}

BOOST_DATA_TEST_CASE_F(any_connection_fixture, connect_any_error, fns_any)
{
    // Setup
    const network_functions_any& fn = sample;

    // Call the function
    fn.connect(conn, connect_params_builder().ssl(ssl_mode::require).database("bad_db").build())
        .validate_error(
            common_server_errc::er_dbaccess_denied_error,
            "Access denied for user 'integ_user'@'%' to database 'bad_db'"
        );
}

// Set character set
BOOST_DATA_TEST_CASE_F(any_connection_fixture, set_character_set_success, fns_any)
{
    // Setup
    const network_functions_any& fn = sample;
    connect();

    // Call the function
    fn.set_character_set(conn, ascii_charset).validate_no_error();

    // Success
    BOOST_TEST(conn.current_character_set()->name == "ascii");
}

BOOST_DATA_TEST_CASE_F(any_connection_fixture, set_character_set_error, fns_any)
{
    // Setup
    const network_functions_any& fn = sample;
    connect();

    // Call the function
    fn.set_character_set(conn, character_set{"bad_charset", nullptr})
        .validate_error(common_server_errc::er_unknown_character_set, "Unknown character set: 'bad_charset'");

    // Success
    BOOST_TEST(conn.current_character_set()->name == "utf8mb4");
}

// Run pipeline
BOOST_DATA_TEST_CASE_F(any_connection_fixture, run_pipeline_success, fns_any)
{
    // Setup
    const network_functions_any& fn = sample;
    pipeline_request req;
    req.add_set_character_set(ascii_charset).add_execute("SET @myvar = 42").add_execute("SELECT @myvar");
    std::vector<stage_response> res;
    connect();

    // Run the function
    fn.run_pipeline(conn, req, res).validate_no_error();

    // Success
    BOOST_TEST(conn.current_character_set().value() == ascii_charset);
    BOOST_TEST_REQUIRE(res.size() == 3u);
    BOOST_TEST(res.at(0).error() == error_code());
    BOOST_TEST(res.at(0).diag() == diagnostics());
    BOOST_TEST(res.at(1).as_results().rows() == rows(), per_element());
    BOOST_TEST(res.at(2).as_results().rows() == makerows(1, 42), per_element());
}

BOOST_DATA_TEST_CASE_F(any_connection_fixture, run_pipeline_error, fns_any)
{
    // Setup
    const network_functions_any& fn = sample;
    pipeline_request req;
    req.add_execute("SET @myvar = 42")
        .add_prepare_statement("SELECT * FROM bad_table")
        .add_execute("SELECT @myvar");
    std::vector<stage_response> res;
    connect();

    // Run the function
    fn.run_pipeline(conn, req, res)
        .validate_error(
            common_server_errc::er_no_such_table,
            "Table 'boost_mysql_integtests.bad_table' doesn't exist"
        );

    // Stages 0 and 2 were executed successfully
    BOOST_TEST(res.size() == 3u);
    BOOST_TEST(res[0].as_results().rows() == rows(), per_element());
    BOOST_TEST(res[1].error() == common_server_errc::er_no_such_table);
    BOOST_TEST(res[1].diag() == create_server_diag("Table 'boost_mysql_integtests.bad_table' doesn't exist"));
    BOOST_TEST(res[2].as_results().rows() == makerows(1, 42), per_element());
}

BOOST_AUTO_TEST_SUITE_END()  // test_spotchecks

}  // namespace