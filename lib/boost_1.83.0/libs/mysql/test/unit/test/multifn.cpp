//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Since integration tests can't reliably test multifunction operations
// that span over multiple messages, we test the complete multifn fllow in this unit tests.

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/connection.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/check_meta.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_multifn)

using test_connection = connection<test_stream>;

using start_query_netm = netfun_maker_mem<void, test_connection, const string_view&, execution_state&>;
using read_resultset_head_netm = netfun_maker_mem<void, test_connection, execution_state&>;
using read_some_rows_netm = netfun_maker_mem<rows_view, test_connection, execution_state&>;

struct
{
    start_query_netm::signature start_execution;
    read_resultset_head_netm::signature read_resultset_head;
    read_some_rows_netm::signature read_some_rows;
    const char* name;
} all_fns[] = {
    {start_query_netm::sync_errc(&test_connection::start_execution),
     read_resultset_head_netm::sync_errc(&test_connection::read_resultset_head),
     read_some_rows_netm::sync_errc(&test_connection::read_some_rows),
     "sync" },
    {start_query_netm::async_errinfo(&test_connection::async_start_execution),
     read_resultset_head_netm::async_errinfo(&test_connection::async_read_resultset_head),
     read_some_rows_netm::async_errinfo(&test_connection::async_read_some_rows),
     "async"},
};

BOOST_AUTO_TEST_CASE(separate_batches)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st;
            test_connection conn;
            conn.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_break()
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
                .add_break()
                .add_bytes(create_text_row_message(3, "abc"))
                .add_break()
                .add_bytes(create_eof_frame(
                    4,
                    ok_builder().affected_rows(10u).info("1st").more_results(true).build()
                ))
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
            fns.start_execution(conn, "SELECT 1", st).validate_no_error();
            BOOST_TEST_REQUIRE(st.should_read_rows());
            check_meta(st.meta(), {column_type::varchar});

            // 1st resultset, row
            auto rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(st.should_read_rows());
            BOOST_TEST(rv == makerows(1, "abc"));

            // 1st resultset, eof
            rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(st.should_read_head());
            BOOST_TEST(rv == makerows(1));
            BOOST_TEST(st.affected_rows() == 10u);
            BOOST_TEST(st.info() == "1st");

            // 2nd resultset, head
            fns.read_resultset_head(conn, st).validate_no_error();
            BOOST_TEST_REQUIRE(st.should_read_rows());
            check_meta(st.meta(), {column_type::decimal});

            // 2nd resultset, row batch
            rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(st.should_read_rows());
            BOOST_TEST(rv == makerows(1, "ab", "plo"));

            // 2nd resultset, last row & eof
            rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(st.complete());
            BOOST_TEST(rv == makerows(1, "hju"));
            BOOST_TEST(st.affected_rows() == 30u);
            BOOST_TEST(st.info() == "2nd");
        }
    }
}

// The server sent us a single, big message with everything
BOOST_AUTO_TEST_CASE(single_read)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st;
            test_connection conn(buffer_params(4096));

            conn.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
                .add_bytes(create_text_row_message(3, "abc"))
                .add_bytes(create_eof_frame(
                    4,
                    ok_builder().affected_rows(10u).info("1st").more_results(true).build()
                ))
                .add_bytes(create_frame(5, {0x01}))
                .add_bytes(create_coldef_frame(6, meta_builder().type(column_type::decimal).build_coldef()))
                .add_bytes(create_text_row_message(7, "ab"))
                .add_bytes(create_text_row_message(8, "plo"))
                .add_bytes(create_text_row_message(9, "hju"))
                .add_bytes(create_eof_frame(10, ok_builder().affected_rows(30u).info("2nd").build()));

            // Start
            fns.start_execution(conn, "SELECT 1", st).validate_no_error();
            BOOST_TEST_REQUIRE(st.should_read_rows());
            check_meta(st.meta(), {column_type::varchar});

            // First resultset
            auto rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(st.should_read_head());
            BOOST_TEST(rv == makerows(1, "abc"));
            BOOST_TEST(st.affected_rows() == 10u);
            BOOST_TEST(st.info() == "1st");

            // 2nd resultset, head
            fns.read_resultset_head(conn, st).validate_no_error();
            BOOST_TEST_REQUIRE(st.should_read_rows());
            check_meta(st.meta(), {column_type::decimal});

            // 2nd resultset
            rv = fns.read_some_rows(conn, st).get();
            BOOST_TEST_REQUIRE(st.complete());
            BOOST_TEST(rv == makerows(1, "ab", "plo", "hju"));
            BOOST_TEST(st.affected_rows() == 30u);
            BOOST_TEST(st.info() == "2nd");
        }
    }
}

BOOST_AUTO_TEST_CASE(empty_resultsets)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            execution_state st;
            test_connection conn(buffer_params(4096));

            conn.stream()
                .add_bytes(
                    create_ok_frame(1, ok_builder().affected_rows(10u).info("1st").more_results(true).build())
                )
                .add_bytes(
                    create_ok_frame(2, ok_builder().affected_rows(20u).info("2nd").more_results(true).build())
                )
                .add_bytes(create_ok_frame(3, ok_builder().affected_rows(30u).info("3rd").build()));

            // Start
            fns.start_execution(conn, "SELECT 1", st).validate_no_error();
            BOOST_TEST_REQUIRE(st.should_read_head());
            BOOST_TEST(st.meta().size() == 0u);
            BOOST_TEST(st.affected_rows() == 10u);
            BOOST_TEST(st.info() == "1st");

            // 2nd resultset
            fns.read_resultset_head(conn, st).validate_no_error();
            BOOST_TEST_REQUIRE(st.should_read_head());
            BOOST_TEST(st.meta().size() == 0u);
            BOOST_TEST(st.affected_rows() == 20u);
            BOOST_TEST(st.info() == "2nd");

            // 3rd resultset
            fns.read_resultset_head(conn, st).validate_no_error();
            BOOST_TEST_REQUIRE(st.complete());
            BOOST_TEST(st.meta().size() == 0u);
            BOOST_TEST(st.affected_rows() == 30u);
            BOOST_TEST(st.info() == "3rd");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
