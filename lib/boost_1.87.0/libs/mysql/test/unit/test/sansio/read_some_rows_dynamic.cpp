//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/mysql/detail/execution_processor/execution_state_impl.hpp>

#include <boost/mysql/impl/internal/sansio/connection_state_data.hpp>
#include <boost/mysql/impl/internal/sansio/read_some_rows_dynamic.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/buffer_concat.hpp"
#include "test_unit/algo_test.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::mysql::detail::execution_state_impl;

BOOST_AUTO_TEST_SUITE(test_read_some_rows_dynamic)

struct fixture : algo_fixture_base
{
    execution_state_impl exec_st;
    detail::read_some_rows_dynamic_algo algo{diag, {&exec_st}};

    fixture()
    {
        // Prepare the state, such that it's ready to read rows
        add_meta(exec_st, {meta_builder().type(column_type::varchar).build_coldef()});
        exec_st.sequence_number() = 42;

        // Put something in shared_fields, simulating a previous read
        st.shared_fields.push_back(field_view("prev"));
    }

    rows_view result() const { return algo.result(st); }
};

BOOST_AUTO_TEST_CASE(eof)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(create_eof_frame(42, ok_builder().affected_rows(1).info("1st").build()))
        .check(fix);

    BOOST_TEST(fix.result() == makerows(1));
    BOOST_TEST_REQUIRE(fix.exec_st.is_complete());
    BOOST_TEST(fix.exec_st.get_affected_rows() == 1u);
    BOOST_TEST(fix.exec_st.get_info() == "1st");
}

BOOST_AUTO_TEST_CASE(batch_with_rows)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(buffer_builder()
                         .add(create_text_row_message(42, "abc"))
                         .add(create_text_row_message(43, "von"))
                         .build())
        .check(fix);

    // Check
    BOOST_TEST(fix.result() == makerows(1, "abc", "von"));
    BOOST_TEST(fix.exec_st.is_reading_rows());
}

BOOST_AUTO_TEST_CASE(batch_with_rows_eof)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(buffer_builder()
                         .add(create_text_row_message(42, "abc"))
                         .add(create_text_row_message(43, "von"))
                         .add(create_eof_frame(44, ok_builder().affected_rows(1).info("1st").build()))
                         .build())
        .check(fix);

    // Check
    BOOST_TEST(fix.result() == makerows(1, "abc", "von"));
    BOOST_TEST_REQUIRE(fix.exec_st.is_complete());
    BOOST_TEST(fix.exec_st.get_affected_rows() == 1u);
    BOOST_TEST(fix.exec_st.get_info() == "1st");
}

// All the other error cases are already tested in read_some_rows_impl. Spotcheck
BOOST_AUTO_TEST_CASE(error)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test().expect_read(client_errc::incomplete_message).check(fix, client_errc::incomplete_message);
}

BOOST_AUTO_TEST_SUITE_END()
