//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/config.hpp>

#ifdef BOOST_MYSQL_CXX14

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connection.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/static_execution_state.hpp>
#include <boost/mysql/static_results.hpp>

#include <boost/core/span.hpp>
#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>
#include <boost/optional/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <tuple>

#include "test_common/check_meta.hpp"
#include "test_common/printing.hpp"
#include "test_integration/common.hpp"
#include "test_integration/metadata_validator.hpp"
#include "test_integration/static_rows.hpp"
#include "test_integration/tcp_network_fixture.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

// Note: the dynamic interface is already covered by stored_procedures, multi_queries, prepared_statements and
// spotchecks

namespace {

BOOST_AUTO_TEST_SUITE(test_static_iface)

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

void validate_multified_rows(boost::span<const row_multifield> r)
{
    BOOST_TEST_REQUIRE(r.size() == 2u);
    BOOST_TEST((r[0] == row_multifield{1.1f, 11, "aaa"}));
    BOOST_TEST((r[1] == row_multifield{{}, 22, "bbb"}));
}

constexpr const char* multifield_bad_msg =
    "NULL checks failed for field 'field_nullable': the database type may be NULL, but the C++ type "
    "cannot. Use std::optional<T> or boost::optional<T>\n"
    "Incompatible types for field 'field_int': C++ type 'string' is not compatible with DB type 'INT'\n"
    "Field 'field_missing' is not present in the data returned by the server";

BOOST_AUTO_TEST_SUITE(singlefn)

BOOST_FIXTURE_TEST_CASE(describe_structs, tcp_network_fixture)
{
    connect();

    static_results<row_multifield> result;
    conn.execute("SELECT * FROM multifield_table ORDER BY id", result);

    // Verify results
    validate_multifield_meta(result.meta());
    validate_multified_rows(result.rows());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(tuples, tcp_network_fixture)
{
    connect();

    using tuple_t = std::tuple<int, std::string, int, boost::optional<float>>;  // trailing fields discarded
    static_results<tuple_t> result;
    conn.execute("SELECT * FROM multifield_table ORDER BY id", result);

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

// This spotchecks having a repeated empty row type, too
BOOST_FIXTURE_TEST_CASE(multi_resultset, tcp_network_fixture)
{
    params.set_multi_queries(true);
    connect();
    start_transaction();

    static_results<row_multifield, empty, row_2fields, empty> result;
    conn.execute(
        R"%(
            SELECT * FROM multifield_table;
            DELETE FROM updates_table;
            SELECT * FROM one_row_table;
            SET @v1 = 2
        )%",
        result
    );

    // Validate results
    validate_multifield_meta(result.meta<0>());
    validate_multified_rows(result.rows<0>());
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

    validate_2fields_meta(result.meta<2>(), "one_row_table");
    BOOST_TEST_REQUIRE(result.rows<2>().size() == 1u);
    BOOST_TEST_REQUIRE((result.rows<2>()[0] == row_2fields{1, std::string("f0")}));
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

BOOST_FIXTURE_TEST_CASE(metadata_check_failed, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;
    static_results<row_multifield_bad> result;
    conn.execute("SELECT * FROM multifield_table ORDER BY id", result, ec, diag);

    BOOST_TEST(ec == client_errc::metadata_check_failed);
    BOOST_TEST(diag.client_message() == multifield_bad_msg);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_empty_resultset, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;
    static_results<std::tuple<int>> result;
    conn.execute("SET @v1 = 2", result, ec, diag);

    const char* expected_msg =
        "Field in position 0 can't be mapped: there are more fields in your C++ data type than in your query";

    BOOST_TEST(ec == client_errc::metadata_check_failed);
    BOOST_TEST(diag.client_message() == expected_msg);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_subsequent_resultset, tcp_network_fixture)
{
    params.set_multi_queries(true);
    connect();

    error_code ec;
    diagnostics diag;
    static_results<empty, row_multifield_bad> result;
    conn.execute("SET @v1 = 2; SELECT * FROM multifield_table ORDER BY id", result, ec, diag);

    BOOST_TEST(ec == client_errc::metadata_check_failed);
    BOOST_TEST(diag.client_message() == multifield_bad_msg);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;
    static_results<row_2fields, empty> result;
    conn.execute("SELECT * FROM one_row_table", result, ec, diag);

    BOOST_TEST(ec == client_errc::num_resultsets_mismatch);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch_empty_resultset, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;
    static_results<empty, empty> result;
    conn.execute("SET @v1 = 2", result, ec, diag);

    BOOST_TEST(ec == client_errc::num_resultsets_mismatch);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(multifn)
BOOST_FIXTURE_TEST_CASE(describe_structs, tcp_network_fixture)
{
    connect();

    // Start
    static_execution_state<row_multifield> result;
    conn.start_execution("SELECT * FROM multifield_table WHERE id = 1", result);
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.should_read_rows());

    // Read rows
    std::array<row_multifield, 3> rws;
    std::size_t num_rows = conn.read_some_rows(result, boost::span<row_multifield>(rws));
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws[0] == row_multifield{1.1f, 11, "aaa"}));

    // Read again, in case the EOF came separately
    num_rows = conn.read_some_rows(result, boost::span<row_multifield>(rws));
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.complete());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(tuples, tcp_network_fixture)
{
    connect();

    using tuple_t = std::tuple<int, std::string, int>;  // trailing fields discarded

    // Start
    static_execution_state<tuple_t> result;
    conn.start_execution("SELECT * FROM multifield_table WHERE id = 1", result);
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.should_read_rows());

    // Read rows
    std::array<tuple_t, 3> rws;
    std::size_t num_rows = conn.read_some_rows(result, boost::span<tuple_t>(rws));
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws[0] == tuple_t{1, "aaa", 11}));

    // Read again, in case the EOF came separately
    num_rows = conn.read_some_rows(result, boost::span<tuple_t>(rws));
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.complete());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

// This spotchecks having repeated empty row types, too
BOOST_FIXTURE_TEST_CASE(multi_resultset, tcp_network_fixture)
{
    params.set_multi_queries(true);
    connect();
    start_transaction();

    // Start
    static_execution_state<row_multifield, empty, row_2fields, empty> result;
    conn.start_execution(
        R"%(
            SELECT * FROM multifield_table WHERE id = 1;
            DELETE FROM updates_table;
            SELECT * FROM one_row_table;
            SET @v1 = 2
        )%",
        result
    );
    validate_multifield_meta(result.meta());
    BOOST_TEST(result.should_read_rows());

    // Read rows (r1)
    std::array<row_multifield, 3> rws;
    std::size_t num_rows = conn.read_some_rows(result, boost::span<row_multifield>(rws));
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws[0] == row_multifield{1.1f, 11, "aaa"}));

    // Read again, in case the EOF came separately (r1)
    num_rows = conn.read_some_rows(result, boost::span<row_multifield>(rws));
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.should_read_head());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Next resultset (r2, empty)
    conn.read_resultset_head(result);
    BOOST_TEST(result.should_read_head());
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 3u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Next resultset (r3)
    conn.read_resultset_head(result);
    BOOST_TEST(result.should_read_rows());
    validate_2fields_meta(result.meta(), "one_row_table");

    // Read rows (r3)
    std::array<row_2fields, 3> rws2;
    num_rows = conn.read_some_rows(result, boost::span<row_2fields>(rws2));
    BOOST_TEST_REQUIRE(num_rows == 1u);
    BOOST_TEST((rws2[0] == row_2fields{1, std::string("f0")}));

    // Read again, in case the EOF came separately (r3)
    num_rows = conn.read_some_rows(result, boost::span<row_2fields>(rws2));
    BOOST_TEST_REQUIRE(num_rows == 0u);
    BOOST_TEST(result.should_read_head());
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");

    // Next resultset (r4, empty)
    conn.read_resultset_head(result);
    BOOST_TEST(result.complete());
    BOOST_TEST(result.meta().size() == 0u);
    BOOST_TEST(result.affected_rows() == 0u);
    BOOST_TEST(result.warning_count() == 0u);
    BOOST_TEST(result.last_insert_id() == 0u);
    BOOST_TEST(result.info() == "");
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;
    static_execution_state<row_multifield_bad> result;
    conn.start_execution("SELECT * FROM multifield_table ORDER BY id", result, ec, diag);

    BOOST_TEST(ec == client_errc::metadata_check_failed);
    BOOST_TEST(diag.client_message() == multifield_bad_msg);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_empty_resultset, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;
    static_execution_state<std::tuple<int>> result;
    conn.start_execution("SET @v1 = 2", result, ec, diag);

    const char* expected_msg =
        "Field in position 0 can't be mapped: there are more fields in your C++ data type than in your query";

    BOOST_TEST(ec == client_errc::metadata_check_failed);
    BOOST_TEST(diag.client_message() == expected_msg);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;
    static_execution_state<row_2fields, empty> result;

    // Start execution
    conn.start_execution("SELECT * FROM empty_table", result);

    // Error is detected when reading the OK packet in read_some_rows
    std::array<row_2fields, 3> storage;
    conn.read_some_rows(result, boost::span<row_2fields>(storage), ec, diag);
    BOOST_TEST(ec == client_errc::num_resultsets_mismatch);
}

BOOST_FIXTURE_TEST_CASE(num_resultsets_mismatch_empty_resultset, tcp_network_fixture)
{
    connect();

    error_code ec;
    diagnostics diag;

    // Start
    static_execution_state<empty, empty> result;
    conn.start_execution("SET @v1 = 2", result, ec, diag);

    BOOST_TEST(ec == client_errc::num_resultsets_mismatch);
}

BOOST_FIXTURE_TEST_CASE(metadata_check_failed_subsequent_resultset, tcp_network_fixture)
{
    params.set_multi_queries(true);
    connect();

    error_code ec;
    diagnostics diag;
    static_execution_state<empty, row_multifield_bad> result;

    // Start execution goes OK
    conn.start_execution("SET @v1 = 2; SELECT * FROM multifield_table", result);

    // Error is detected when reading next head
    conn.read_resultset_head(result, ec, diag);

    BOOST_TEST(ec == client_errc::metadata_check_failed);
    BOOST_TEST(diag.client_message() == multifield_bad_msg);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace

#endif
