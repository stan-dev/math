//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connection.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/mysql/detail/config.hpp>

#include "test_common/create_basic.hpp"
#include "test_integration/common.hpp"
#include "test_integration/er_connection.hpp"
#include "test_integration/metadata_validator.hpp"
#include "test_integration/static_rows.hpp"

using namespace boost::mysql::test;
using boost::mysql::common_server_errc;
using boost::mysql::execution_state;
using boost::mysql::field_view;
using boost::mysql::results;
using boost::mysql::row_view;

namespace {

BOOST_AUTO_TEST_SUITE(test_spotchecks)

auto err_net_samples = create_network_samples({
    "tcp_sync_errc",
    "tcp_sync_exc",
    "tcp_async_callback",
    "tcp_async_coroutines",
});

// Handshake
BOOST_MYSQL_NETWORK_TEST(handshake_success, network_fixture, all_network_samples())
{
    setup_and_physical_connect(sample.net);
    conn->handshake(params).validate_no_error();
    BOOST_TEST(conn->uses_ssl() == var->supports_ssl());
}

BOOST_MYSQL_NETWORK_TEST(handshake_error, network_fixture, err_net_samples)
{
    setup_and_physical_connect(sample.net);
    params.set_database("bad_database");
    conn->handshake(params).validate_error(
        common_server_errc::er_dbaccess_denied_error,
        {"database", "bad_database"}
    );
}

// Connect: success is already widely tested throughout integ tests
BOOST_MYSQL_NETWORK_TEST(connect_error, network_fixture, err_net_samples)
{
    setup(sample.net);
    set_credentials("integ_user", "bad_password");
    conn->connect(params).validate_error(
        common_server_errc::er_access_denied_error,
        {"access denied", "integ_user"}
    );
    BOOST_TEST(!conn->is_open());
}

// Start query (legacy)
BOOST_MYSQL_NETWORK_TEST(start_query_legacy_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    execution_state st;
    conn->start_query("SELECT * FROM empty_table", st).get();
    BOOST_TEST(st.should_read_rows());
    validate_2fields_meta(st.meta(), "empty_table");
}

BOOST_MYSQL_NETWORK_TEST(start_query_legacy_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);

    execution_state st;
    conn->start_query("SELECT field_varchar, field_bad FROM one_row_table", st)
        .validate_error(common_server_errc::er_bad_field_error, {"unknown column", "field_bad"});
}

// Start execution (query)
BOOST_MYSQL_NETWORK_TEST(start_execution_query_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    execution_state st;
    conn->start_execution("SELECT * FROM empty_table", st).get();
    BOOST_TEST(st.should_read_rows());
    validate_2fields_meta(st.meta(), "empty_table");
}

BOOST_MYSQL_NETWORK_TEST(start_execution_query_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);

    execution_state st;
    conn->start_execution("SELECT field_varchar, field_bad FROM one_row_table", st)
        .validate_error(common_server_errc::er_bad_field_error, {"unknown column", "field_bad"});
}

// Query (legacy)
BOOST_MYSQL_NETWORK_TEST(query_legacy_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    results result;
    conn->query("SELECT 'hello', 42", result).get();
    BOOST_TEST(result.rows().size() == 1u);
    BOOST_TEST(result.rows()[0] == makerow("hello", 42));
    BOOST_TEST(result.meta().size() == 2u);
}

BOOST_MYSQL_NETWORK_TEST(query_legacy_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);

    results result;
    conn->query("SELECT field_varchar, field_bad FROM one_row_table", result)
        .validate_error(common_server_errc::er_bad_field_error, {"unknown column", "field_bad"});
}

// execute (query)
BOOST_MYSQL_NETWORK_TEST(execute_query_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    results result;
    conn->execute("SELECT 'hello', 42", result).get();
    BOOST_TEST(result.rows().size() == 1u);
    BOOST_TEST(result.rows()[0] == makerow("hello", 42));
    BOOST_TEST(result.meta().size() == 2u);
}

BOOST_MYSQL_NETWORK_TEST(execute_query_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);

    results result;
    conn->execute("SELECT field_varchar, field_bad FROM one_row_table", result)
        .validate_error(common_server_errc::er_bad_field_error, {"unknown column", "field_bad"});
}

// Prepare statement
BOOST_MYSQL_NETWORK_TEST(prepare_statement_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();
    BOOST_TEST_REQUIRE(stmt.valid());
    BOOST_TEST(stmt.id() > 0u);
    BOOST_TEST(stmt.num_params() == 2u);
}

BOOST_MYSQL_NETWORK_TEST(prepare_statement_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);
    conn->prepare_statement("SELECT * FROM bad_table WHERE id IN (?, ?)")
        .validate_error(common_server_errc::er_no_such_table, {"table", "doesn't exist", "bad_table"});
}

// Start statement execution (legacy, iterator)
BOOST_MYSQL_NETWORK_TEST(start_statement_execution_legacy_it_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Prepare
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Execute
    execution_state st;
    std::forward_list<field_view> params{field_view("item"), field_view(42)};
    conn->start_statement_execution(stmt, params.begin(), params.end(), st).validate_no_error();
    validate_2fields_meta(st.meta(), "empty_table");
    BOOST_TEST(st.should_read_rows());
}

BOOST_MYSQL_NETWORK_TEST(start_statement_execution_legacy_it_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);
    start_transaction();

    // Prepare
    auto stmt = conn->prepare_statement("INSERT INTO inserts_table (field_varchar, field_date) VALUES (?, ?)")
                    .get();

    // Execute
    execution_state st;
    std::forward_list<field_view> params{field_view("f0"), field_view("bad_date")};
    conn->start_statement_execution(stmt, params.begin(), params.end(), st)
        .validate_error(
            common_server_errc::er_truncated_wrong_value,
            {"field_date", "bad_date", "incorrect date value"}
        );
}

// Start execution (statement, iterator). No error spotcheck, since it's the same underlying function
BOOST_MYSQL_NETWORK_TEST(start_execution_stmt_it_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Prepare
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Execute
    execution_state st;
    std::forward_list<field_view> params{field_view("item"), field_view(42)};
    conn->start_execution(stmt.bind(params.cbegin(), params.cend()), st).validate_no_error();
    validate_2fields_meta(st.meta(), "empty_table");
    BOOST_TEST(st.should_read_rows());
}

// Start statement execution (legacy, tuple)
BOOST_MYSQL_NETWORK_TEST(
    start_statement_execution_legacy_tuple_success,
    network_fixture,
    all_network_samples()
)
{
    setup_and_connect(sample.net);

    // Prepare
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Execute
    execution_state st;
    conn->start_statement_execution(stmt, field_view(42), field_view(40), st).validate_no_error();
    validate_2fields_meta(st.meta(), "empty_table");
    BOOST_TEST(st.should_read_rows());
}

BOOST_MYSQL_NETWORK_TEST(start_statement_execution_legacy_tuple_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);
    start_transaction();

    // Prepare
    auto stmt = conn->prepare_statement("INSERT INTO inserts_table (field_varchar, field_date) VALUES (?, ?)")
                    .get();

    // Execute
    execution_state st;
    conn->start_statement_execution(stmt, field_view("abc"), field_view("bad_date"), st)
        .validate_error(
            common_server_errc::er_truncated_wrong_value,
            {"field_date", "bad_date", "incorrect date value"}
        );
}

// start execution (statement, tuple). No error spotcheck since it's the same underlying fn
BOOST_MYSQL_NETWORK_TEST(start_execution_statement_tuple_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Prepare
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Execute
    execution_state st;
    conn->start_execution(stmt.bind(field_view(42), field_view(40)), st).validate_no_error();
    validate_2fields_meta(st.meta(), "empty_table");
    BOOST_TEST(st.should_read_rows());
}

// Execute statement (legacy)
BOOST_MYSQL_NETWORK_TEST(execute_statement_legacy_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Prepare
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Execute
    results result;
    conn->execute_statement(stmt, field_view("item"), field_view(42), result).validate_no_error();
    BOOST_TEST(result.rows().size() == 0u);
}

BOOST_MYSQL_NETWORK_TEST(execute_statement_legacy_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);
    start_transaction();

    // Prepare
    auto stmt = conn->prepare_statement("INSERT INTO inserts_table (field_varchar, field_date) VALUES (?, ?)")
                    .get();

    // Execute
    results result;
    conn->execute_statement(stmt, field_view("f0"), field_view("bad_date"), result)
        .validate_error(
            common_server_errc::er_truncated_wrong_value,
            {"field_date", "bad_date", "incorrect date value"}
        );
}

// Execute (statement, iterator). No error spotcheck since it's the same underlying fn
BOOST_MYSQL_NETWORK_TEST(execute_statement_iterator_success, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);

    // Prepare
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Execute
    results result;
    std::forward_list<field_view> params{field_view("item"), field_view(42)};
    conn->execute(stmt.bind(params.cbegin(), params.cend()), result).validate_no_error();
    BOOST_TEST(result.rows().size() == 0u);
}

// Execute (statement, tuple). No error spotcheck since it's the same underlying fn
BOOST_MYSQL_NETWORK_TEST(execute_statement_tuple_success, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);

    // Prepare
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Execute
    results result;
    conn->execute(stmt.bind(field_view("item"), field_view(42)), result).validate_no_error();
    BOOST_TEST(result.rows().size() == 0u);
}

// Close statement: no server error spotcheck
BOOST_MYSQL_NETWORK_TEST(close_statement_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Prepare a statement
    auto stmt = conn->prepare_statement("SELECT * FROM empty_table WHERE id IN (?, ?)").get();

    // Close the statement
    conn->close_statement(stmt).validate_no_error();

    // The statement is no longer valid
    results result;
    conn->execute_statement(stmt, field_view("a"), field_view("b"), result).validate_any_error();
}

// Read some rows: no server error spotcheck
BOOST_MYSQL_NETWORK_TEST(read_some_rows_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Generate an execution state
    execution_state st;
    conn->start_query("SELECT * FROM one_row_table", st);
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // Read once. st may or may not be complete, depending
    // on how the buffer reallocated memory
    auto rows = conn->read_some_rows(st).get();
    BOOST_TEST((rows == makerows(2, 1, "f0")));

    // Reading again should complete st
    rows = conn->read_some_rows(st).get();
    BOOST_TEST(rows.empty());
    validate_eof(st);

    // Reading again does nothing
    rows = conn->read_some_rows(st).get();
    BOOST_TEST(rows.empty());
    validate_eof(st);
}

// Read resultset head
BOOST_MYSQL_NETWORK_TEST(read_resultset_head_success, network_fixture, all_network_samples())
{
    params.set_multi_queries(true);
    setup_and_connect(sample.net);

    // Generate an execution state
    execution_state st;
    conn->start_query("SELECT * FROM empty_table; SELECT * FROM one_row_table", st);
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // Read the OK packet to finish 1st resultset
    conn->read_some_rows(st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_head());

    // Read head
    conn->read_resultset_head(st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // Reading head again does nothing
    conn->read_resultset_head(st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // We can read rows now
    auto rows = conn->read_some_rows(st).get();
    BOOST_TEST((rows == makerows(2, 1, "f0")));
}

BOOST_MYSQL_NETWORK_TEST(read_resultset_head_error, network_fixture, all_network_samples())
{
    params.set_multi_queries(true);
    setup_and_connect(sample.net);

    // Generate an execution state
    execution_state st;
    conn->start_query("SELECT * FROM empty_table; SELECT bad_field FROM one_row_table", st);
    BOOST_TEST_REQUIRE(st.should_read_rows());

    // Read the OK packet to finish 1st resultset
    conn->read_some_rows(st).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_head());

    // Read head for the 2nd resultset. This one contains an error, which is detected when reading head.
    conn->read_resultset_head(st).validate_error_exact(
        common_server_errc::er_bad_field_error,
        "Unknown column 'bad_field' in 'field list'"
    );
}

// Ping
BOOST_MYSQL_NETWORK_TEST(ping_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);
    conn->ping().validate_no_error();
}

BOOST_MYSQL_NETWORK_TEST(ping_error, network_fixture, all_network_samples())
{
    setup(sample.net);

    // Ping should return an error for an unconnected connection
    conn->ping().validate_any_error();
}

// Quit connection: no server error spotcheck
BOOST_MYSQL_NETWORK_TEST(quit_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Quit
    conn->quit().validate_no_error();

    // We are no longer able to query
    results result;
    conn->query("SELECT 1", result).validate_any_error();
}

// Close connection: no server error spotcheck
BOOST_MYSQL_NETWORK_TEST(close_connection_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    // Close
    conn->close().validate_no_error();

    // We are no longer able to query
    boost::mysql::results result;
    conn->query("SELECT 1", result).validate_any_error();

    // The stream is closed
    BOOST_TEST(!conn->is_open());

    // Closing again returns OK (and does nothing)
    conn->close().validate_no_error();

    // Stream (socket) still closed
    BOOST_TEST(!conn->is_open());
}

// TODO: move this to a unit test
BOOST_MYSQL_NETWORK_TEST(not_open_connection, network_fixture, err_net_samples)
{
    setup(sample.net);
    conn->close().validate_no_error();
    BOOST_TEST(!conn->is_open());
}

#ifdef BOOST_MYSQL_CXX14
// Execute (static) - errors are already covered by the other tests
BOOST_MYSQL_NETWORK_TEST(execute_static_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    er_connection::static_results_t result;
    conn->execute("CALL sp_spotchecks()", result).get();
    BOOST_TEST(result.rows<0>().size() == 1u);
    BOOST_TEST((result.rows<0>()[0] == row_multifield{1.1f, 11, "aaa"}));
}

// start_execution, read_resultset_head, read_some_rows success
BOOST_MYSQL_NETWORK_TEST(start_execution_and_followups_static_success, network_fixture, all_network_samples())
{
    setup_and_connect(sample.net);

    er_connection::static_state_t st;

    // Start
    conn->start_execution("CALL sp_spotchecks()", st).get();
    BOOST_TEST(st.should_read_rows());

    // Read r1 rows
    std::array<row_multifield, 2> storage;
    std::size_t num_rows = conn->read_some_rows(st, storage).get();
    BOOST_TEST(num_rows == 1u);
    BOOST_TEST((storage[0] == row_multifield{1.1f, 11, "aaa"}));

    // Ensure we're in the next resultset
    num_rows = conn->read_some_rows(st, storage).get();
    BOOST_TEST(num_rows == 0u);
    BOOST_TEST(st.should_read_head());

    // Read r2 head
    conn->read_resultset_head(st).get();
    BOOST_TEST(st.should_read_rows());

    // Read r2 rows
    std::array<row_2fields, 2> storage2;
    num_rows = conn->read_some_rows(st, storage2).get();
    BOOST_TEST(num_rows == 1u);
    BOOST_TEST((storage2[0] == row_2fields{1, std::string("f0")}));

    // Ensure we're in the next resultset
    num_rows = conn->read_some_rows(st, storage2).get();
    BOOST_TEST(num_rows == 0u);
    BOOST_TEST(st.should_read_head());

    // Read r3 head (empty)
    conn->read_resultset_head(st).get();
    BOOST_TEST(st.complete());
}

// read_some_rows failure (the other error cases are already widely tested)
BOOST_MYSQL_NETWORK_TEST(read_some_rows_error, network_fixture, err_net_samples)
{
    setup_and_connect(sample.net);

    er_connection::static_state_t st;

    // Start
    conn->start_execution("SELECT * FROM multifield_table WHERE id = 42", st).get();
    BOOST_TEST(st.should_read_rows());

    // No rows matched, so reading rows reads the OK packet. This will report the num resultsets mismatch
    std::array<row_multifield, 2> storage;
    conn->read_some_rows(st, storage)
        .validate_error_exact_client(boost::mysql::client_errc::num_resultsets_mismatch);
}
#endif

BOOST_AUTO_TEST_SUITE_END()  // test_spotchecks

}  // namespace
