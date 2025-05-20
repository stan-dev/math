//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/results.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/check_meta.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/metadata_validator.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::test_tools::per_element;

namespace {

BOOST_AUTO_TEST_SUITE(test_stored_procedures)

BOOST_FIXTURE_TEST_CASE(without_selects, any_connection_fixture)
{
    connect();
    start_transaction();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_insert(?)", as_netresult).get();

    // Call the procedure
    results result;
    conn.async_execute(stmt.bind("abc"), result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST_REQUIRE(result.size() == 1u);
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.rows() == rows(), per_element());
    BOOST_TEST(result.affected_rows() == 1u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);  // this refers to the CALL, not to the INSERT
    BOOST_TEST(result.info() == "");
    BOOST_TEST(result.out_params() == row_view());

    // Verify it took place
    conn.async_execute("SELECT field_varchar FROM inserts_table", result, as_netresult).validate_no_error();
    BOOST_TEST(result.rows() == makerows(1, "abc"), per_element());
}

BOOST_FIXTURE_TEST_CASE(with_one_select, any_connection_fixture)
{
    connect();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_select_1(?)", as_netresult).get();

    // Call the procedure
    results result;
    conn.async_execute(stmt.bind("abc"), result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST_REQUIRE(result.size() == 2u);
    validate_2fields_meta(result[0].meta(), "one_row_table");
    BOOST_TEST(result[0].rows() == makerows(2, 1, "f0"), per_element());
    BOOST_TEST(result[0].affected_rows() == 0u);
    BOOST_TEST(result[0].warning_count() == 0u);
    BOOST_TEST(result[0].last_insert_id() == 0u);
    BOOST_TEST(result[0].info() == "");
    BOOST_TEST(result[1].meta().size() == 0u);
    BOOST_TEST(result[1].rows() == rows(), per_element());
    BOOST_TEST(result[1].affected_rows() == 0u);
    BOOST_TEST(result[1].warning_count() == 0u);
    BOOST_TEST(result[1].last_insert_id() == 0u);
    BOOST_TEST(result[1].info() == "");
    BOOST_TEST(result.out_params() == row_view());
}

BOOST_FIXTURE_TEST_CASE(with_two_selects, any_connection_fixture)
{
    connect();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_select_2(?, ?)", as_netresult).get();

    // Call the procedure
    results result;
    conn.async_execute(stmt.bind("abc", 42), result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST_REQUIRE(result.size() == 3u);
    validate_2fields_meta(result[0].meta(), "one_row_table");
    BOOST_TEST(result[0].rows() == makerows(2, 1, "f0"), per_element());
    BOOST_TEST(result[0].affected_rows() == 0u);
    BOOST_TEST(result[0].warning_count() == 0u);
    BOOST_TEST(result[0].last_insert_id() == 0u);
    BOOST_TEST(result[0].info() == "");
    check_meta(result[1].meta(), {column_type::varchar, column_type::int_});
    BOOST_TEST(result[1].rows() == makerows(2, "abc", 42), per_element());
    BOOST_TEST(result[1].affected_rows() == 0u);
    BOOST_TEST(result[1].warning_count() == 0u);
    BOOST_TEST(result[1].last_insert_id() == 0u);
    BOOST_TEST(result[1].info() == "");
    BOOST_TEST(result[2].meta().size() == 0u);
    BOOST_TEST(result[2].rows() == rows(), per_element());
    BOOST_TEST(result[2].affected_rows() == 0u);
    BOOST_TEST(result[2].warning_count() == 0u);
    BOOST_TEST(result[2].last_insert_id() == 0u);
    BOOST_TEST(result[2].info() == "");
    BOOST_TEST(result.out_params() == row_view());
}

BOOST_FIXTURE_TEST_CASE(with_two_selects_multifn, any_connection_fixture)
{
    connect();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_select_2(?, ?)", as_netresult).get();

    // Call the procedure
    execution_state st;
    conn.async_start_execution(stmt.bind("abc", 42), st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    validate_2fields_meta(st.meta(), "one_row_table");

    // Read rows for the 1st select
    auto rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(rv == makerows(2, 1, "f0"), per_element());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
    BOOST_TEST(!st.is_out_params());

    // Read 2nd resultset's head
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    check_meta(st.meta(), {column_type::varchar, column_type::int_});

    // Read 2nd resultset's rows
    rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(rv == makerows(2, "abc", 42), per_element());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
    BOOST_TEST(!st.is_out_params());

    // Read final resultset
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.meta().size() == 0u);
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
    BOOST_TEST(!st.is_out_params());
}

BOOST_FIXTURE_TEST_CASE(output_params_not_bound, any_connection_fixture)
{
    connect();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_outparams(?, @var1, @var2)", as_netresult).get();

    // Call the procedure
    results result;
    conn.async_execute(stmt.bind(10), result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST_REQUIRE(result.size() == 2u);
    validate_2fields_meta(result[0].meta(), "one_row_table");
    BOOST_TEST(result[0].rows() == makerows(2, 1, "f0"), per_element());
    BOOST_TEST(result[1].meta().size() == 0u);
    BOOST_TEST(result[1].rows() == rows(), per_element());
    BOOST_TEST(result[1].affected_rows() == 0u);
    BOOST_TEST(result[1].warning_count() == 0u);
    BOOST_TEST(result[1].last_insert_id() == 0u);
    BOOST_TEST(result[1].info() == "");
    BOOST_TEST(result.out_params() == row_view());
}

BOOST_FIXTURE_TEST_CASE(output_params_bound, any_connection_fixture)
{
    connect();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_outparams(?, ?, ?)", as_netresult).get();

    // Call the procedure
    results result;
    conn.async_execute(stmt.bind(10, nullptr, 30), result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST_REQUIRE(result.size() == 3u);
    validate_2fields_meta(result[0].meta(), "one_row_table");
    BOOST_TEST(result[0].rows() == makerows(2, 1, "f0"), per_element());
    BOOST_TEST(!result[0].is_out_params());
    check_meta(result[1].meta(), {column_type::int_, column_type::int_});
    BOOST_TEST(result[1].rows() == makerows(2, 10, 31), per_element());
    BOOST_TEST(result[1].affected_rows() == 0u);
    BOOST_TEST(result[1].warning_count() == 0u);
    BOOST_TEST(result[1].last_insert_id() == 0u);
    BOOST_TEST(result[1].info() == "");
    BOOST_TEST(result[1].is_out_params());
    BOOST_TEST(result[2].meta().size() == 0u);
    BOOST_TEST(result[2].rows() == rows(), per_element());
    BOOST_TEST(result[2].affected_rows() == 0u);
    BOOST_TEST(result[2].warning_count() == 0u);
    BOOST_TEST(result[2].last_insert_id() == 0u);
    BOOST_TEST(result[2].info() == "");
    BOOST_TEST(!result[2].is_out_params());
    BOOST_TEST(result.out_params() == makerow(10, 31));
}

BOOST_FIXTURE_TEST_CASE(output_params_bound_multifn, any_connection_fixture)
{
    connect();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_outparams(?, ?, ?)", as_netresult).get();

    // Call the procedure
    execution_state st;
    conn.async_start_execution(stmt.bind(10, nullptr, 30), st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    validate_2fields_meta(st.meta(), "one_row_table");

    // 1st resultset, rows
    auto rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST(rv == makerows(2, 1, "f0"), per_element());
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(!st.is_out_params());

    // out params, head
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    check_meta(st.meta(), {column_type::int_, column_type::int_});

    // out params, rows and eof
    rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(rv == makerows(2, 10, 31), per_element());
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
    BOOST_TEST(st.is_out_params());

    // final eof
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.meta().size() == 0u);
    BOOST_TEST(st.affected_rows() == 0u);
    BOOST_TEST(st.warning_count() == 0u);
    BOOST_TEST(st.last_insert_id() == 0u);
    BOOST_TEST(st.info() == "");
    BOOST_TEST(!st.is_out_params());
}

BOOST_FIXTURE_TEST_CASE(with_signal, any_connection_fixture)
{
    connect();

    // Statement
    auto stmt = conn.async_prepare_statement("CALL sp_signal()", as_netresult).get();

    // Call the procedure. It should fail, since we're invoking SIGNAL
    results result;
    conn.async_execute(stmt.bind(), result, as_netresult)
        .validate_error(common_server_errc::er_no, "An error occurred");
}

BOOST_FIXTURE_TEST_CASE(with_query, any_connection_fixture)
{
    connect();

    // Call the procedure
    results result;
    conn.async_execute("CALL sp_outparams(42, @var1, @var2)", result, as_netresult).validate_no_error();

    // Verify results
    BOOST_TEST_REQUIRE(result.size() == 2u);
    validate_2fields_meta(result[0].meta(), "one_row_table");
    BOOST_TEST(result[0].rows() == makerows(2, 1, "f0"));
    BOOST_TEST(result[1].meta().size() == 0u);
    BOOST_TEST(result[1].rows() == rows(), per_element());
    BOOST_TEST(result[1].affected_rows() == 0u);
    BOOST_TEST(result[1].warning_count() == 0u);
    BOOST_TEST(result[1].last_insert_id() == 0u);
    BOOST_TEST(result[1].info() == "");
    BOOST_TEST(result.out_params() == row_view());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
