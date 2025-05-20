//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/execution_processor/execution_processor.hpp>
#include <boost/mysql/detail/execution_processor/execution_state_impl.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include "execution_processor_helpers.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/printing.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::mysql::detail::execution_state_impl;
using boost::mysql::detail::output_ref;
using boost::mysql::detail::resultset_encoding;

namespace {

BOOST_AUTO_TEST_SUITE(test_execution_state_impl)

// OK packet checking
void check_ok_r1(const execution_state_impl& st)
{
    BOOST_TEST(st.get_affected_rows() == 1u);
    BOOST_TEST(st.get_last_insert_id() == 2u);
    BOOST_TEST(st.get_warning_count() == 4u);
    BOOST_TEST(st.get_info() == "Information");
    BOOST_TEST(st.get_is_out_params() == false);
}

void check_ok_r2(const execution_state_impl& st)
{
    BOOST_TEST(st.get_affected_rows() == 5u);
    BOOST_TEST(st.get_last_insert_id() == 6u);
    BOOST_TEST(st.get_warning_count() == 8u);
    BOOST_TEST(st.get_info() == "more_info");
    BOOST_TEST(st.get_is_out_params() == true);
}

void check_ok_r3(const execution_state_impl& st)
{
    BOOST_TEST(st.get_affected_rows() == 10u);
    BOOST_TEST(st.get_last_insert_id() == 11u);
    BOOST_TEST(st.get_warning_count() == 12u);
    BOOST_TEST(st.get_info() == "");
    BOOST_TEST(st.get_is_out_params() == false);
}

struct fixture
{
    std::vector<field_view> fields;
    execution_state_impl st;
    boost::mysql::diagnostics diag;
};

BOOST_FIXTURE_TEST_CASE(one_resultset_data, fixture)
{
    // Initial. Verify that we clear any previous result
    exec_access(st)
        .reset(resultset_encoding::binary)
        .meta({column_type::geometry})
        .ok(ok_builder().affected_rows(1).last_insert_id(2).warnings(3).more_results(true).info("abc").build()
        );

    // Reset
    st.reset(resultset_encoding::text, metadata_mode::full);
    BOOST_TEST(st.is_reading_first());

    // Head indicates resultset with metadata
    st.on_num_meta(2);
    BOOST_TEST(st.is_reading_meta());

    // First metadata
    auto err = st.on_meta(create_meta_r1_0(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_meta());

    // Second metadata, ready to read rows
    err = st.on_meta(create_meta_r1_1(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    check_meta_r1(st.meta());

    // Rows
    auto r1 = create_text_row_body(10, "abc");
    auto r2 = create_text_row_body(20, "cdef");
    err = st.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    err = st.on_row(r2, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(fields == make_fv_vector(10, "abc", 20, "cdef"));

    // End of resultset
    err = st.on_row_ok_packet(create_ok_r1());
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_r1(st.meta());
    check_ok_r1(st);
}

BOOST_FIXTURE_TEST_CASE(one_resultset_empty, fixture)
{
    // Directly end of resultet, no meta
    auto err = st.on_head_ok_packet(create_ok_r1(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_empty(st.meta());
    check_ok_r1(st);
}

BOOST_FIXTURE_TEST_CASE(two_resultsets_data_data, fixture)
{
    // Resultset r1 (rows are not stored anyhow in execution states)
    add_meta(st, create_meta_r1());

    // OK packet indicates more results
    auto err = st.on_row_ok_packet(create_ok_r1(true));
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_first_subseq());
    check_meta_r1(st.meta());
    check_ok_r1(st);

    // Resultset r2: indicates resultset with meta
    st.on_num_meta(1);
    BOOST_TEST(st.is_reading_meta());

    // First packet
    err = st.on_meta(create_meta_r2_0(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    check_meta_r2(st.meta());

    // Rows
    auto r1 = create_text_row_body(90u);
    err = st.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    BOOST_TEST(fields == make_fv_vector(90u));

    // OK packet, no more resultsets
    err = st.on_row_ok_packet(create_ok_r2());
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_r2(st.meta());
    check_ok_r2(st);
}

BOOST_FIXTURE_TEST_CASE(two_resultsets_empty_data, fixture)
{
    // Resultset r1
    add_ok(st, create_ok_r1(true));
    BOOST_TEST(st.is_reading_first_subseq());
    check_meta_empty(st.meta());
    check_ok_r1(st);

    // Resultset r2: indicates data
    st.on_num_meta(1);
    BOOST_TEST(st.is_reading_meta());

    // Metadata packet
    auto err = st.on_meta(create_meta_r2_0(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    check_meta_r2(st.meta());

    // Rows
    auto r1 = create_text_row_body(90u);
    auto r2 = create_text_row_body(100u);
    err = st.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    err = st.on_row(r2, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());

    // Final OK packet
    err = st.on_row_ok_packet(create_ok_r2());
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_r2(st.meta());
    check_ok_r2(st);
}

BOOST_FIXTURE_TEST_CASE(two_resultsets_data_empty, fixture)
{
    // Resultset r1
    add_meta(st, create_meta_r1());

    // OK packet indicates more results
    auto err = st.on_row_ok_packet(create_ok_r1(true));
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_first_subseq());
    check_meta_r1(st.meta());
    check_ok_r1(st);

    // OK packet for 2nd result
    err = st.on_head_ok_packet(create_ok_r2(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_empty(st.meta());
    check_ok_r2(st);
}

BOOST_FIXTURE_TEST_CASE(two_resultsets_empty_empty, fixture)
{
    // OK packet indicates more results
    auto err = st.on_head_ok_packet(create_ok_r1(true), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_first_subseq());
    check_meta_empty(st.meta());
    check_ok_r1(st);

    // OK packet for 2nd result
    err = st.on_head_ok_packet(create_ok_r2(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_empty(st.meta());
    check_ok_r2(st);
}

BOOST_FIXTURE_TEST_CASE(three_resultsets_empty_empty_data, fixture)
{
    // Two first resultsets
    add_ok(st, create_ok_r1(true));
    auto err = st.on_head_ok_packet(create_ok_r2(true), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_first_subseq());
    check_meta_empty(st.meta());
    check_ok_r2(st);

    // Resultset r3: head indicates resultset with metadata
    st.on_num_meta(3);
    BOOST_TEST(st.is_reading_meta());

    // Metadata
    err = st.on_meta(create_meta_r3_0(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_meta());

    err = st.on_meta(create_meta_r3_1(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_meta());

    err = st.on_meta(create_meta_r3_2(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    check_meta_r3(st.meta());

    // Rows
    auto r1 = create_text_row_body(4.2f, 90.0, 9);
    err = st.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    BOOST_TEST(fields == make_fv_vector(4.2f, 90.0, 9));

    // End of resultset
    err = st.on_row_ok_packet(create_ok_r3());
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_r3(st.meta());
    check_ok_r3(st);
}

BOOST_FIXTURE_TEST_CASE(three_resultsets_data_empty_data, fixture)
{
    // Two first resultsets
    exec_access(st).meta(create_meta_r1()).ok(create_ok_r1(true));
    auto err = st.on_head_ok_packet(create_ok_r2(true), diag);
    BOOST_TEST(st.is_reading_first_subseq());
    check_meta_empty(st.meta());
    check_ok_r2(st);

    // Resultset r3: head indicates resultset with metadata
    st.on_num_meta(3);
    BOOST_TEST(st.is_reading_meta());

    // Metadata
    err = st.on_meta(create_meta_r3_0(), diag);
    throw_on_error(err, diag);
    err = st.on_meta(create_meta_r3_1(), diag);
    throw_on_error(err, diag);
    err = st.on_meta(create_meta_r3_2(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(st.is_reading_rows());
    check_meta_r3(st.meta());

    // Rows
    auto r1 = create_text_row_body(4.2f, 90.0, 9);
    err = st.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(fields == make_fv_vector(4.2f, 90.0, 9));

    // End of resultset
    err = st.on_row_ok_packet(create_ok_r3());
    throw_on_error(err, diag);
    BOOST_TEST(st.is_complete());
    check_meta_r3(st.meta());
    check_ok_r3(st);
}

BOOST_FIXTURE_TEST_CASE(info_string_ownserhip, fixture)
{
    // OK packet received, doesn't own the string
    std::string info = "Some info";
    auto err = st.on_head_ok_packet(ok_builder().more_results(true).info(info).build(), diag);
    throw_on_error(err, diag);

    // st does, so changing info doesn't affect
    info = "other info";
    BOOST_TEST(st.get_info() == "Some info");

    // Repeat the process for row OK packet
    st.on_num_meta(1);
    err = st.on_meta(meta_builder().build_coldef(), diag);
    throw_on_error(err, diag);
    err = st.on_row_ok_packet(ok_builder().info(info).build());
    throw_on_error(err, diag);
    info = "abcdfefgh";
    BOOST_TEST(st.get_info() == "other info");
}

BOOST_FIXTURE_TEST_CASE(error_deserializing_row, fixture)
{
    add_meta(st, create_meta_r1());
    auto bad_row = create_text_row_body(42, "abc");
    bad_row.push_back(0xff);

    auto err = st.on_row(bad_row, output_ref(), fields);

    BOOST_TEST(err == client_errc::extra_bytes);
}

BOOST_FIXTURE_TEST_CASE(meta_mode_minimal, fixture)
{
    st.reset(resultset_encoding::text, metadata_mode::minimal);
    add_meta(st, create_meta_r1());
    BOOST_TEST(st.meta()[0].column_name() == "");
}

BOOST_FIXTURE_TEST_CASE(meta_mode_full, fixture)
{
    st.reset(resultset_encoding::text, metadata_mode::full);
    add_meta(st, create_meta_r1());
    BOOST_TEST(st.meta()[0].column_name() == "ftiny");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace