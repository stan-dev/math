//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Since integration tests can't reliably test multifunction operations
// that span over multiple messages, we test the complete multifn fllow in this unit tests.

#include <boost/mysql/any_connection.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/check_meta.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/network_result.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/test_any_connection.hpp"
#include "test_unit/test_stream.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_multifn)

BOOST_FIXTURE_TEST_CASE(separate_batches, io_context_fixture)
{
    execution_state st;
    auto conn = create_test_any_connection(ctx);
    get_stream(conn)
        .add_bytes(create_frame(1, {0x01}))
        .add_break()
        .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .add_break()
        .add_bytes(create_text_row_message(3, "abc"))
        .add_break()
        .add_bytes(create_eof_frame(4, ok_builder().affected_rows(10u).info("1st").more_results(true).build())
        )
        .add_break()
        .add_bytes(create_frame(5, {0x01}))
        .add_break()
        .add_bytes(create_coldef_frame(6, meta_builder().type(column_type::decimal).build_coldef()))
        .add_break()
        .add_bytes(create_text_row_message(7, "ab"))
        .add_bytes(create_text_row_message(8, "plo"))
        .add_break()
        .add_bytes(create_text_row_message(9, "hju"))
        .add_bytes(create_eof_frame(10, ok_builder().affected_rows(30u).info("2nd").build()));

    // Start
    conn.async_start_execution("SELECT 1", st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    check_meta(st.meta(), {column_type::varchar});

    // 1st resultset, row
    auto rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    BOOST_TEST(rv == makerows(1, "abc"));

    // 1st resultset, eof
    rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(rv == makerows(1));
    BOOST_TEST(st.affected_rows() == 10u);
    BOOST_TEST(st.info() == "1st");

    // 2nd resultset, head
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    check_meta(st.meta(), {column_type::decimal});

    // 2nd resultset, row batch
    rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    BOOST_TEST(rv == makerows(1, "ab", "plo"));

    // 2nd resultset, last row & eof
    rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(rv == makerows(1, "hju"));
    BOOST_TEST(st.affected_rows() == 30u);
    BOOST_TEST(st.info() == "2nd");
}

// The server sent us a single, big message with everything
BOOST_FIXTURE_TEST_CASE(single_read, io_context_fixture)
{
    execution_state st;
    any_connection_params params;
    params.initial_buffer_size = 4096;
    auto conn = create_test_any_connection(ctx, params);

    get_stream(conn)
        .add_bytes(create_frame(1, {0x01}))
        .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .add_bytes(create_text_row_message(3, "abc"))
        .add_bytes(create_eof_frame(4, ok_builder().affected_rows(10u).info("1st").more_results(true).build())
        )
        .add_bytes(create_frame(5, {0x01}))
        .add_bytes(create_coldef_frame(6, meta_builder().type(column_type::decimal).build_coldef()))
        .add_bytes(create_text_row_message(7, "ab"))
        .add_bytes(create_text_row_message(8, "plo"))
        .add_bytes(create_text_row_message(9, "hju"))
        .add_bytes(create_eof_frame(10, ok_builder().affected_rows(30u).info("2nd").build()));

    // Start
    conn.async_start_execution("SELECT 1", st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    check_meta(st.meta(), {column_type::varchar});

    // First resultset
    auto rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(rv == makerows(1, "abc"));
    BOOST_TEST(st.affected_rows() == 10u);
    BOOST_TEST(st.info() == "1st");

    // 2nd resultset, head
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_rows());
    check_meta(st.meta(), {column_type::decimal});

    // 2nd resultset
    rv = conn.async_read_some_rows(st, as_netresult).get();
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(rv == makerows(1, "ab", "plo", "hju"));
    BOOST_TEST(st.affected_rows() == 30u);
    BOOST_TEST(st.info() == "2nd");
}

BOOST_FIXTURE_TEST_CASE(empty_resultsets, io_context_fixture)
{
    execution_state st;
    any_connection_params params;
    params.initial_buffer_size = 4096;
    auto conn = create_test_any_connection(ctx, params);

    get_stream(conn)
        .add_bytes(create_ok_frame(1, ok_builder().affected_rows(10u).info("1st").more_results(true).build()))
        .add_bytes(create_ok_frame(2, ok_builder().affected_rows(20u).info("2nd").more_results(true).build()))
        .add_bytes(create_ok_frame(3, ok_builder().affected_rows(30u).info("3rd").build()));

    // Start
    conn.async_start_execution("SELECT 1", st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(st.meta().size() == 0u);
    BOOST_TEST(st.affected_rows() == 10u);
    BOOST_TEST(st.info() == "1st");

    // 2nd resultset
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.should_read_head());
    BOOST_TEST(st.meta().size() == 0u);
    BOOST_TEST(st.affected_rows() == 20u);
    BOOST_TEST(st.info() == "2nd");

    // 3rd resultset
    conn.async_read_resultset_head(st, as_netresult).validate_no_error();
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.meta().size() == 0u);
    BOOST_TEST(st.affected_rows() == 30u);
    BOOST_TEST(st.info() == "3rd");
}

BOOST_AUTO_TEST_SUITE_END()
