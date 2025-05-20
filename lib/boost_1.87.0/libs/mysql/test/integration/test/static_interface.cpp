//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/config.hpp>

#ifdef BOOST_MYSQL_CXX14
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/pfr.hpp>
#include <boost/mysql/static_execution_state.hpp>
#include <boost/mysql/static_results.hpp>

#include <boost/core/span.hpp>
#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>
#include <boost/optional/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/pfr/ops.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstdint>
#include <tuple>

#include "test_common/check_meta.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/metadata_validator.hpp"
#include "test_integration/static_rows.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::span;
using boost::test_tools::per_element;

// Note: the dynamic interface is already covered by stored_procedures, multi_queries, prepared_statements and
// spotchecks

BOOST_AUTO_TEST_SUITE(test_static_iface)

// A row type like row_multifield, but without Describe metadata
struct row_multifield_pfr
{
    boost::optional<float> field_nullable;
    std::int32_t field_int;
    std::string field_varchar;
};

// Same, but only with literal fields (required by PFR in C++14 mode)
struct row_multifield_pfr_literal
{
    std::int64_t id;
    std::int32_t field_int;
    double field_double;
};

void validate_multifield_meta(metadata_collection_view meta)
{
    check_meta(
        meta,
        {
            column_type::int_,
            column_type::varchar,
            column_type::int_,
            column_type::float_,
            column_type::double_,
        }
    );
}

std::array<row_multifield, 2> expected_multifield_rows()
{
    return {
        {
         {1.1f, 11, "aaa"},
         {{}, 22, "bbb"},
         }
    };
}

constexpr const char* multifield_bad_msg =
    "NULL checks failed for field 'field_nullable': the database type may be NULL, but the C++ type "
    "cannot. Use std::optional<T> or boost::optional<T>\n"
    "Incompatible types for field 'field_int': C++ type 'string' is not compatible with DB type 'INT'\n"
    "Field 'field_missing' is not present in the data returned by the server";

BOOST_AUTO_TEST_SUITE(singlefn)

BOOST_FIXTURE_TEST_CASE(describe_structs, any_connection_fixture)
{
    connect();

    static_results<row_multifield> result;
    conn.async_execute("SELECT * FROM multifield_table ORDER BY id", result, as_netresult)
        .validate_no_error();

    // Verify results
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.rows() == expected_multifield_rows(), per_element());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(tuples, any_connection_fixture)
{
    connect();

    using tuple_t = std::tuple<int, std::string, int, boost::optional<float>>;  // trailing fields discarded
    static_results<tuple_t> result;
    conn.async_execute("SELECT * FROM multifield_table ORDER BY id", result, as_netresult)
        .validate_no_error();

    // Verify results
    validate_multifield_meta(result.meta());
    BOOST_TEST_REQUIRE(result.rows().size() == 2u);
    BOOST_TEST((result.rows()[0] == tuple_t{1, "aaa", 11, 1.1f}));
    BOOST_TEST((result.rows()[1] == tuple_t{2, "bbb", 22, {}}));
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

#if BOOST_PFR_CORE_NAME_ENABLED
BOOST_FIXTURE_TEST_CASE(pfr_structs_by_name, any_connection_fixture)
{
    connect();

    static_results<pfr_by_name<row_multifield_pfr>> result;
    conn.async_execute("SELECT * FROM multifield_table ORDER BY id", result, as_netresult)
        .validate_no_error();

    // Verify results
    validate_multifield_meta(result.meta());
    BOOST_TEST_REQUIRE(result.rows().size() == 2u);
    BOOST_TEST(boost::pfr::eq(result.rows()[0], row_multifield_pfr{1.1f, 11, "aaa"}));
    BOOST_TEST(boost::pfr::eq(result.rows()[1], row_multifield_pfr{{}, 22, "bbb"}));
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}
#endif

#if BOOST_PFR_ENABLED && defined(BOOST_MYSQL_CXX14)
BOOST_FIXTURE_TEST_CASE(pfr_structs_by_position, any_connection_fixture)
{
    connect();

    static_results<pfr_by_position<row_multifield_pfr_literal>> result;
    conn.async_execute(
            "SELECT id, field_int, field_double FROM multifield_table ORDER BY id",
            result,
            as_netresult
    )
        .validate_no_error();

    // Verify results
    check_meta(result.meta(), {column_type::int_, column_type::int_, column_type::double_});
    BOOST_TEST_REQUIRE(result.rows().size() == 2u);
    BOOST_TEST(boost::pfr::eq(result.rows()[0], row_multifield_pfr_literal{1, 11, 0.1}));
    BOOST_TEST(boost::pfr::eq(result.rows()[1], row_multifield_pfr_literal{2, 22, 0.2}));
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}
#endif

// This spotchecks having a repeated empty row type, too
BOOST_FIXTURE_TEST_CASE(multi_resultset, any_connection_fixture)
{
    connect(connect_params_builder().multi_queries(true).build());
    start_transaction();

    static_results<row_multifield, empty, row_2fields, empty> result;
    constexpr const char* query =
        "SELECT * FROM multifield_table;"
        "DELETE FROM updates_table;"
        "SELECT * FROM one_row_table;"
        "SET @v1 = 2";
    conn.async_execute(query, result, as_netresult).validate_no_error();

    // Validate results
    validate_multifield_meta(result.meta<0>());
    BOOST_TEST(result.rows<0>() == expected_multifield_rows(), per_element());
    BOOST_TEST(result.affected_rows<0>() == 0u);
    BOOST_TEST(result.warning_count<0>() == 0u);
    BOOST_TEST(result.last_insert_id<0>() == 0u);
    BOOST_TEST(result.info<0>() == "");

    BOOST_TEST(result.meta<1>().size() == 0u);
    BOOST_TEST(result.rows<1>().size() == 0u);
    BOOST_TEST(result.affected_rows<1>() == 3u);
    BOOST_TEST(result.warning_count<1>() == 0u);
    BOOST_TEST(result.last_insert_id<1>() == 0u);
    BOOST_TEST(result.info<1>() == "");

    const row_2fields expected_2fields[]{
        {1, std::string("f0")}
    };
    validate_2fields_meta(result.meta<2>(), "one_row_table");
    BOOST_TEST(result.rows<2>() == expected_2fields, per_element());
    BOOST_TEST(result.affected_rows<2>() == 0u);
    BOOST_TEST(result.warning_count<2>() == 0u);
    BOOST_TEST(result.last_insert_id<2>() == 0u);
    BOOST_TEST(result.info<2>() == "");

    BOOST_TEST(result.meta<3>().size() == 0u);
    BOOST_TEST(result.rows<3>().size() == 0u);
    BOOST_TEST(result.affected_rows<3>() == 0u);
    BOOST_TEST(result.warning_count<3>() == 0u);
    BOOST_TEST(result.last_insert_id<3>() == 0u);
    BOOST_TEST(result.info<3>() == "");
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed, any_connection_fixture)
{
    connect();

    static_results<row_multifield_bad> result;
    conn.async_execute("SELECT * FROM multifield_table ORDER BY id", result, as_netresult)
        .validate_error(client_errc::metadata_check_failed, multifield_bad_msg);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_empty_resultset, any_connection_fixture)
{
    connect();

    const char* expected_msg =
        "Field in position 0 can't be mapped: there are more fields in your C++ data type than in your query";
    static_results<std::tuple<int>> result;
    conn.async_execute("SET @v1 = 2", result, as_netresult)
        .validate_error(client_errc::metadata_check_failed, expected_msg);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_subsequent_resultset, any_connection_fixture)
{
    connect(connect_params_builder().multi_queries(true).build());

    static_results<empty, row_multifield_bad> result;
    conn.async_execute("SET @v1 = 2; SELECT * FROM multifield_table ORDER BY id", result, as_netresult)
        .validate_error(client_errc::metadata_check_failed, multifield_bad_msg);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch, any_connection_fixture)
{
    connect();

    static_results<row_2fields, empty> result;
    conn.async_execute("SELECT * FROM one_row_table", result, as_netresult)
        .validate_error(client_errc::num_resultsets_mismatch);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch_empty_resultset, any_connection_fixture)
{
    connect();

    static_results<empty, empty> result;
    conn.async_execute("SET @v1 = 2", result, as_netresult)
        .validate_error(client_errc::num_resultsets_mismatch);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(multifn)
BOOST_FIXTURE_TEST_CASE(describe_structs, any_connection_fixture)
{
    connect();

    // Start
    static_execution_state<row_multifield> result;
    conn.async_start_execution("SELECT * FROM multifield_table WHERE id = 1", result, as_netresult)
        .validate_no_error();
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.should_read_rows());

    // Read rows
    std::array<row_multifield, 3> rws;
    std::size_t num_rows = conn.async_read_some_rows(result, span<row_multifield>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws[0] == row_multifield{1.1f, 11, "aaa"}));

    // Read again, in case the EOF came separately
    num_rows = conn.async_read_some_rows(result, span<row_multifield>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.complete());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(tuples, any_connection_fixture)
{
    connect();

    using tuple_t = std::tuple<int, std::string, int>;  // trailing fields discarded

    // Start
    static_execution_state<tuple_t> result;
    conn.async_start_execution("SELECT * FROM multifield_table WHERE id = 1", result, as_netresult)
        .validate_no_error();
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.should_read_rows());

    // Read rows
    std::array<tuple_t, 3> rws;
    std::size_t num_rows = conn.async_read_some_rows(result, span<tuple_t>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws[0] == tuple_t{1, "aaa", 11}));

    // Read again, in case the EOF came separately
    num_rows = conn.async_read_some_rows(result, span<tuple_t>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.complete());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

#if BOOST_PFR_CORE_NAME_ENABLED
BOOST_FIXTURE_TEST_CASE(pfr_structs_by_name, any_connection_fixture)
{
    connect();

    // Start
    static_execution_state<pfr_by_name<row_multifield_pfr>> result;
    conn.async_start_execution("SELECT * FROM multifield_table WHERE id = 1", result, as_netresult)
        .validate_no_error();
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.should_read_rows());

    // Read rows
    std::array<row_multifield_pfr, 3> rws;
    std::size_t num_rows = conn.async_read_some_rows(result, span<row_multifield_pfr>(rws), as_netresult)
                               .get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST(boost::pfr::eq(rws[0], row_multifield_pfr{1.1f, 11, "aaa"}));

    // Read again, in case the EOF came separately
    num_rows = conn.async_read_some_rows(result, span<row_multifield_pfr>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.complete());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}
#endif

#if BOOST_PFR_ENABLED && defined(BOOST_MYSQL_CXX14)
BOOST_FIXTURE_TEST_CASE(pfr_structs_by_position, any_connection_fixture)
{
    connect();

    // Start
    static_execution_state<pfr_by_position<row_multifield_pfr_literal>> result;
    conn.async_start_execution(
            "SELECT id, field_int, field_double FROM multifield_table WHERE id = 1",
            result,
            as_netresult
    )
        .validate_no_error();
    check_meta(result.meta(), {column_type::int_, column_type::int_, column_type::double_});
    BOOST_TEST(result.should_read_rows());

    // Read rows
    std::array<row_multifield_pfr_literal, 3> rws;
    auto num_rows = conn.async_read_some_rows(result, span<row_multifield_pfr_literal>(rws), as_netresult)
                        .get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST(boost::pfr::eq(rws[0], row_multifield_pfr_literal{1, 11, 0.1}));

    // Read again, in case the EOF came separately
    num_rows = conn.async_read_some_rows(result, span<row_multifield_pfr_literal>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.complete());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}
#endif

// This spotchecks having repeated empty row types, too
BOOST_FIXTURE_TEST_CASE(multi_resultset, any_connection_fixture)
{
    connect(connect_params_builder().multi_queries(true).build());
    start_transaction();

    // Start
    constexpr const char* query =
        "SELECT * FROM multifield_table WHERE id = 1;"
        "DELETE FROM updates_table;"
        "SELECT * FROM one_row_table;"
        "SET @v1 = 2";
    static_execution_state<row_multifield, empty, row_2fields, empty> result;
    conn.async_start_execution(query, result, as_netresult).validate_no_error();
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.should_read_rows());

    // Read rows (r1)
    std::array<row_multifield, 3> rws;
    std::size_t num_rows = conn.async_read_some_rows(result, span<row_multifield>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws[0] == row_multifield{1.1f, 11, "aaa"}));

    // Read again, in case the EOF came separately (r1)
    num_rows = conn.async_read_some_rows(result, span<row_multifield>(rws), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.should_read_head());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Next resultset (r2, empty)
    conn.async_read_resultset_head(result, as_netresult).validate_no_error();
    BOOST_TEST(result.should_read_head());
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 3u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Next resultset (r3)
    conn.async_read_resultset_head(result, as_netresult).validate_no_error();
    BOOST_TEST(result.should_read_rows());
    validate_2fields_meta(result.meta(), "one_row_table");

    // Read rows (r3)
    std::array<row_2fields, 3> rws2;
    num_rows = conn.async_read_some_rows(result, span<row_2fields>(rws2), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws2[0] == row_2fields{1, std::string("f0")}));

    // Read again, in case the EOF came separately (r3)
    num_rows = conn.async_read_some_rows(result, span<row_2fields>(rws2), as_netresult).get();
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.should_read_head());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Next resultset (r4, empty)
    conn.async_read_resultset_head(result, as_netresult).validate_no_error();
    BOOST_TEST(result.complete());
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed, any_connection_fixture)
{
    connect();

    static_execution_state<row_multifield_bad> result;
    conn.async_start_execution("SELECT * FROM multifield_table ORDER BY id", result, as_netresult)
        .validate_error(client_errc::metadata_check_failed, multifield_bad_msg);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_empty_resultset, any_connection_fixture)
{
    connect();

    const char* expected_msg =
        "Field in position 0 can't be mapped: there are more fields in your C++ data type than in your query";
    static_execution_state<std::tuple<int>> result;
    conn.async_start_execution("SET @v1 = 2", result, as_netresult)
        .validate_error(client_errc::metadata_check_failed, expected_msg);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch, any_connection_fixture)
{
    connect();

    static_execution_state<row_2fields, empty> result;

    // Start execution
    conn.async_start_execution("SELECT * FROM empty_table", result, as_netresult).validate_no_error();

    // Error is detected when reading the OK packet in read_some_rows
    std::array<row_2fields, 3> storage;
    conn.async_read_some_rows(result, span<row_2fields>(storage), as_netresult)
        .validate_error(client_errc::num_resultsets_mismatch);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch_empty_resultset, any_connection_fixture)
{
    connect();

    // Start
    static_execution_state<empty, empty> result;
    conn.async_start_execution("SET @v1 = 2", result, as_netresult)
        .validate_error(client_errc::num_resultsets_mismatch);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_subsequent_resultset, any_connection_fixture)
{
    connect(connect_params_builder().multi_queries(true).build());

    static_execution_state<empty, row_multifield_bad> result;

    // Start execution goes OK
    conn.async_start_execution("SET @v1 = 2; SELECT * FROM multifield_table", result, as_netresult)
        .validate_no_error();

    // Error is detected when reading next head
    conn.async_read_resultset_head(result, as_netresult)
        .validate_error(client_errc::metadata_check_failed, multifield_bad_msg);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

#endif
