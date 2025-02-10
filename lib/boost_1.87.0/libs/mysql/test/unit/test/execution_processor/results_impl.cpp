//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/throw_on_error.hpp>

#include <boost/mysql/detail/execution_processor/execution_processor.hpp>
#include <boost/mysql/detail/execution_processor/results_impl.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <boost/test/unit_test.hpp>

#include "execution_processor_helpers.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/printing.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::mysql::detail::output_ref;
using boost::mysql::detail::results_impl;
using boost::mysql::detail::resultset_container;
using boost::mysql::detail::resultset_encoding;

namespace {

BOOST_AUTO_TEST_SUITE(test_results_impl)

// OK packet checking
void check_ok_r1(const results_impl& st, std::size_t idx)
{
    BOOST_TEST(st.get_affected_rows(idx) == 1u);
    BOOST_TEST(st.get_last_insert_id(idx) == 2u);
    BOOST_TEST(st.get_warning_count(idx) == 4u);
    BOOST_TEST(st.get_info(idx) == "Information");
    BOOST_TEST(st.get_is_out_params(idx) == false);
}

void check_ok_r2(const results_impl& st, std::size_t idx)
{
    BOOST_TEST(st.get_affected_rows(idx) == 5u);
    BOOST_TEST(st.get_last_insert_id(idx) == 6u);
    BOOST_TEST(st.get_warning_count(idx) == 8u);
    BOOST_TEST(st.get_info(idx) == "more_info");
    BOOST_TEST(st.get_is_out_params(idx) == true);
}

void check_ok_r3(const results_impl& st, std::size_t idx)
{
    BOOST_TEST(st.get_affected_rows(idx) == 10u);
    BOOST_TEST(st.get_last_insert_id(idx) == 11u);
    BOOST_TEST(st.get_warning_count(idx) == 12u);
    BOOST_TEST(st.get_info(idx) == "");
    BOOST_TEST(st.get_is_out_params(idx) == false);
}

BOOST_AUTO_TEST_SUITE(resultset_container_)
BOOST_AUTO_TEST_CASE(append_from_empty)
{
    // Initial
    resultset_container c;
    BOOST_TEST(c.empty());
    BOOST_TEST(c.size() == 0u);

    // Append first
    c.emplace_back().num_rows = 1;
    BOOST_TEST(!c.empty());
    BOOST_TEST(c.size() == 1u);
    BOOST_TEST(c[0].num_rows == 1u);
    BOOST_TEST(c.back().num_rows == 1u);

    // Append second
    c.emplace_back().num_rows = 2;
    BOOST_TEST(!c.empty());
    BOOST_TEST(c.size() == 2u);
    BOOST_TEST(c[0].num_rows == 1u);
    BOOST_TEST(c[1].num_rows == 2u);
    BOOST_TEST(c.back().num_rows == 2u);

    // Append third
    c.emplace_back().num_rows = 3;
    BOOST_TEST(!c.empty());
    BOOST_TEST(c.size() == 3u);
    BOOST_TEST(c[0].num_rows == 1u);
    BOOST_TEST(c[1].num_rows == 2u);
    BOOST_TEST(c[2].num_rows == 3u);
    BOOST_TEST(c.back().num_rows == 3u);
}

BOOST_AUTO_TEST_CASE(append_from_cleared)
{
    // Initial
    resultset_container c;
    c.emplace_back().num_rows = 42u;
    c.emplace_back().num_rows = 43u;

    // Clear
    c.clear();
    BOOST_TEST(c.empty());
    BOOST_TEST(c.size() == 0u);

    // Append first
    c.emplace_back().num_rows = 1;
    BOOST_TEST(!c.empty());
    BOOST_TEST(c.size() == 1u);
    BOOST_TEST(c[0].num_rows == 1u);
    BOOST_TEST(c.back().num_rows == 1u);

    // Append second
    c.emplace_back().num_rows = 2;
    BOOST_TEST(!c.empty());
    BOOST_TEST(c.size() == 2u);
    BOOST_TEST(c[0].num_rows == 1u);
    BOOST_TEST(c[1].num_rows == 2u);
    BOOST_TEST(c.back().num_rows == 2u);

    // Append third
    c.emplace_back().num_rows = 3;
    BOOST_TEST(!c.empty());
    BOOST_TEST(c.size() == 3u);
    BOOST_TEST(c[0].num_rows == 1u);
    BOOST_TEST(c[1].num_rows == 2u);
    BOOST_TEST(c[2].num_rows == 3u);
    BOOST_TEST(c.back().num_rows == 3u);
}

BOOST_AUTO_TEST_CASE(clear_empty)
{
    // Initial
    resultset_container c;
    c.clear();
    BOOST_TEST(c.empty());
    BOOST_TEST(c.size() == 0u);
}

BOOST_AUTO_TEST_CASE(several_clears)
{
    // Initial
    resultset_container c;
    c.emplace_back().num_rows = 42u;

    // Clear
    c.clear();
    BOOST_TEST(c.empty());
    BOOST_TEST(c.size() == 0u);

    // Append again
    c.emplace_back().num_rows = 1;
    c.emplace_back().num_rows = 2;

    // Clear again
    c.clear();
    BOOST_TEST(c.empty());
    BOOST_TEST(c.size() == 0u);
}
BOOST_AUTO_TEST_SUITE_END()

struct fixture
{
    diagnostics diag;
    std::vector<field_view> fields;
    results_impl r;
};

BOOST_FIXTURE_TEST_CASE(one_resultset_data, fixture)
{
    // Initial. Check that we reset any previous state
    exec_access(r)
        .meta({column_type::geometry})
        .row(makebv("\0\0"))
        .row(makebv("abc"))
        .ok(ok_builder().affected_rows(40).info("some_info").more_results(true).build())
        .meta({column_type::varchar, column_type::mediumint})
        .row("aaaa", 42)
        .ok(ok_builder().info("more_info").more_results(true).build());
    r.reset(resultset_encoding::text, metadata_mode::minimal);
    BOOST_TEST(r.is_reading_first());

    // Head indicates resultset with two columns
    r.on_num_meta(2);
    BOOST_TEST(r.is_reading_meta());

    // First meta
    auto err = r.on_meta(create_meta_r1_0(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_meta());

    // Second meta, ready to read rows
    err = r.on_meta(create_meta_r1_1(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_rows());

    // Rows
    auto r1 = create_text_row_body(42, "abc");
    r.on_row_batch_start();
    err = r.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_rows());

    // End of resultset
    err = r.on_row_ok_packet(create_ok_r1());
    throw_on_error(err, diag);
    r.on_row_batch_finish();  // EOF is part of the batch

    // Verify results
    BOOST_TEST(r.is_complete());
    check_meta_r1(r.get_meta(0));
    check_ok_r1(r, 0);
    BOOST_TEST(r.num_resultsets() == 1u);
    BOOST_TEST(r.get_rows(0) == makerows(2, 42, "abc"));
    BOOST_TEST(r.get_out_params() == row_view());
    BOOST_TEST(fields.empty());  // unused
}

BOOST_FIXTURE_TEST_CASE(one_resultset_empty, fixture)
{
    // Initial
    BOOST_TEST(r.is_reading_first());

    // End of resultset
    auto err = r.on_head_ok_packet(create_ok_r1(), diag);
    throw_on_error(err, diag);

    // verify
    BOOST_TEST(r.is_complete());
    check_meta_empty(r.get_meta(0));
    check_ok_r1(r, 0);
    BOOST_TEST(r.num_resultsets() == 1u);
    BOOST_TEST(r.get_rows(0) == rows_view());
    BOOST_TEST(r.get_out_params() == row_view());
}

BOOST_FIXTURE_TEST_CASE(two_resultsets_data_data, fixture)
{
    // Resultset r1
    exec_access(r).meta(create_meta_r1()).row(42, "abc").row(50, "def");

    // OK packet indicates more results
    auto err = r.on_row_ok_packet(create_ok_r1(true));
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_first_subseq());

    // Resultset r2: indicates resultset with meta
    r.on_num_meta(1);
    BOOST_TEST(r.is_reading_meta());

    // Meta
    err = r.on_meta(create_meta_r2_0(), diag);
    BOOST_TEST(r.is_reading_rows());

    // Row
    auto r1 = create_text_row_body(70);
    r.on_row_batch_start();
    err = r.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);

    // OK packet, no more resultsets
    err = r.on_row_ok_packet(create_ok_r2());
    throw_on_error(err, diag);
    r.on_row_batch_finish();

    // Verify
    BOOST_TEST(r.is_complete());
    check_meta_r1(r.get_meta(0));
    check_meta_r2(r.get_meta(1));
    check_ok_r1(r, 0);
    check_ok_r2(r, 1);
    BOOST_TEST(r.num_resultsets() == 2u);
    BOOST_TEST(r.get_rows(0) == makerows(2, 42, "abc", 50, "def"));
    BOOST_TEST(r.get_rows(1) == makerows(1, 70));
    BOOST_TEST(r.get_out_params() == makerow(70));
}

BOOST_FIXTURE_TEST_CASE(two_resultsets_empty_data, fixture)
{
    // Empty resultset r1, indicating more results
    auto err = r.on_head_ok_packet(create_ok_r1(true), diag);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_first_subseq());

    // Resultset r2: indicates data
    r.on_num_meta(1);
    BOOST_TEST(r.is_reading_meta());

    // Metadata packet
    err = r.on_meta(create_meta_r2_0(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_rows());

    // Rows
    auto r1 = create_text_row_body(70);
    r.on_row_batch_start();
    err = r.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_rows());

    // Final OK packet
    err = r.on_row_ok_packet(create_ok_r2());
    throw_on_error(err, diag);
    r.on_row_batch_finish();

    // Verify
    BOOST_TEST(r.is_complete());
    check_meta_empty(r.get_meta(0));
    check_meta_r2(r.get_meta(1));
    check_ok_r1(r, 0);
    check_ok_r2(r, 1);
    BOOST_TEST(r.num_resultsets() == 2u);
    BOOST_TEST(r.get_rows(0) == rows_view());
    BOOST_TEST(r.get_rows(1) == makerows(1, 70));
    BOOST_TEST(r.get_out_params() == makerow(70));
}

// Note: this tests also an edge case where a resultset indicates
// that contains OUT parameters but is empty
BOOST_FIXTURE_TEST_CASE(two_resultsets_data_empty, fixture)
{
    // Resultset r1
    exec_access(r).meta(create_meta_r1()).row(42, "abc").row(50, "def");

    // OK packet indicates more results
    auto err = r.on_row_ok_packet(create_ok_r1(true));
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_first_subseq());

    // OK packet for 2nd result
    err = r.on_head_ok_packet(create_ok_r2(), diag);
    throw_on_error(err, diag);

    // Verify
    BOOST_TEST(r.is_complete());
    check_meta_r1(r.get_meta(0));
    check_meta_empty(r.get_meta(1));
    check_ok_r1(r, 0);
    check_ok_r2(r, 1);
    BOOST_TEST(r.num_resultsets() == 2u);
    BOOST_TEST(r.get_rows(0) == makerows(2, 42, "abc", 50, "def"));
    BOOST_TEST(r.get_rows(1) == rows_view());
    BOOST_TEST(r.get_out_params() == row_view());
}

BOOST_FIXTURE_TEST_CASE(two_resultsets_empty_empty, fixture)
{
    // Resultset r1
    auto err = r.on_head_ok_packet(create_ok_r1(true), diag);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_first_subseq());

    // OK packet for 2nd result
    err = r.on_head_ok_packet(create_ok_r2(), diag);
    throw_on_error(err, diag);

    // Verify
    BOOST_TEST(r.is_complete());
    check_meta_empty(r.get_meta(0));
    check_meta_empty(r.get_meta(1));
    check_ok_r1(r, 0);
    check_ok_r2(r, 1);
    BOOST_TEST(r.num_resultsets() == 2u);
    BOOST_TEST(r.get_rows(0) == rows_view());
    BOOST_TEST(r.get_rows(1) == rows_view());
    BOOST_TEST(r.get_out_params() == row_view());
}

BOOST_FIXTURE_TEST_CASE(three_resultsets_empty_empty_data, fixture)
{
    // First resultset
    add_ok(r, create_ok_r1(true));

    // Second resultset: OK packet indicates more results
    auto err = r.on_head_ok_packet(create_ok_r2(true), diag);
    BOOST_TEST(r.is_reading_first_subseq());

    // Resultset r3: head indicates resultset with metadata
    r.on_num_meta(3);
    BOOST_TEST(r.is_reading_meta());

    // Metadata
    err = r.on_meta(create_meta_r3_0(), diag);
    throw_on_error(err, diag);
    err = r.on_meta(create_meta_r3_1(), diag);
    throw_on_error(err, diag);
    err = r.on_meta(create_meta_r3_2(), diag);
    throw_on_error(err, diag);
    BOOST_TEST(r.is_reading_rows());

    // Read rows
    auto r1 = create_text_row_body(4.2f, 5.0, 8);
    auto r2 = create_text_row_body(42.0f, 50.0, 80);
    r.on_row_batch_start();
    err = r.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    err = r.on_row(r2, output_ref(), fields);
    throw_on_error(err, diag);

    // End of resultset
    err = r.on_row_ok_packet(create_ok_r3());
    throw_on_error(err, diag);
    r.on_row_batch_finish();

    // Verify
    BOOST_TEST(r.is_complete());
    check_meta_empty(r.get_meta(0));
    check_meta_empty(r.get_meta(1));
    check_meta_r3(r.get_meta(2));
    check_ok_r1(r, 0);
    check_ok_r2(r, 1);
    check_ok_r3(r, 2);
    BOOST_TEST(r.num_resultsets() == 3u);
    BOOST_TEST(r.get_rows(0) == rows_view());
    BOOST_TEST(r.get_rows(1) == rows_view());
    BOOST_TEST(r.get_rows(2) == makerows(3, 4.2f, 5.0, 8, 42.0f, 50.0, 80));
    BOOST_TEST(r.get_out_params() == row_view());
}

// Verify that we do row slicing correctly
BOOST_FIXTURE_TEST_CASE(three_resultsets_data_data_data, fixture)
{
    // Two first resultets
    exec_access(r)
        .meta(create_meta_r1())
        .row(42, "abc")
        .row(50, "def")
        .ok(create_ok_r1(true))
        .meta(create_meta_r2())
        .row(60);

    // OK packet indicates more results
    auto err = r.on_row_ok_packet(create_ok_r2(true));
    throw_on_error(err, diag);

    // Third resultset meta
    r.on_num_meta(3);
    err = r.on_meta(create_meta_r3_0(), diag);
    throw_on_error(err, diag);
    err = r.on_meta(create_meta_r3_1(), diag);
    throw_on_error(err, diag);
    err = r.on_meta(create_meta_r3_2(), diag);
    throw_on_error(err, diag);

    // Rows
    auto r1 = create_text_row_body(4.2f, 5.0, 8);
    auto r2 = create_text_row_body(42.0f, 50.0, 80);
    r.on_row_batch_start();
    err = r.on_row(r1, output_ref(), fields);
    throw_on_error(err, diag);
    err = r.on_row(r2, output_ref(), fields);
    throw_on_error(err, diag);
    r.on_row_batch_finish();

    // OK packet
    err = r.on_row_ok_packet(create_ok_r3());
    throw_on_error(err, diag);

    // Check results
    BOOST_TEST(r.is_complete());
    check_meta_r1(r.get_meta(0));
    check_meta_r2(r.get_meta(1));
    check_meta_r3(r.get_meta(2));
    check_ok_r1(r, 0);
    check_ok_r2(r, 1);
    check_ok_r3(r, 2);
    BOOST_TEST(r.num_resultsets() == 3u);
    BOOST_TEST(r.get_rows(0) == makerows(2, 42, "abc", 50, "def"));
    BOOST_TEST(r.get_rows(1) == makerows(1, 60));
    BOOST_TEST(r.get_rows(2) == makerows(3, 4.2f, 5.0, 8, 42.0f, 50.0, 80));
    BOOST_TEST(r.get_out_params() == makerow(60));
}

BOOST_FIXTURE_TEST_CASE(info_string_ownserhip, fixture)
{
    // Head OK packet
    std::string info = "Some info";
    auto err = r.on_head_ok_packet(ok_builder().more_results(true).info(info).build(), diag);
    throw_on_error(err, diag);

    // Empty OK packet
    info = "";
    err = r.on_head_ok_packet(ok_builder().more_results(true).info(info).build(), diag);

    // Row OK packet
    info = "other info";
    add_meta(r, create_meta_r2());
    err = r.on_row_ok_packet(ok_builder().info(info).build());
    info = "abcdfefgh";
    BOOST_TEST(r.get_info(0) == "Some info");
    BOOST_TEST(r.get_info(1) == "");
    BOOST_TEST(r.get_info(2) == "other info");
}

BOOST_FIXTURE_TEST_CASE(multiple_row_batches, fixture)
{
    // Initial
    add_meta(r, create_meta_r1());

    // Buffers
    auto r1 = create_text_row_body(42, "abc");
    auto r2 = create_text_row_body(50, "bdef");
    auto r3 = create_text_row_body(60, "pov");

    // First batch
    r.on_row_batch_start();
    auto err = r.on_row(r1, output_ref(), fields);
    throw_on_error(err);
    err = r.on_row(r2, output_ref(), fields);
    throw_on_error(err);
    r.on_row_batch_finish();

    // Second batch (only one row)
    r.on_row_batch_start();
    err = r.on_row(r3, output_ref(), fields);
    throw_on_error(err);

    // End of resultset
    err = r.on_row_ok_packet(create_ok_r1());
    throw_on_error(err);
    r.on_row_batch_finish();

    // Verify
    BOOST_TEST(r.is_complete());
    BOOST_TEST(r.num_resultsets() == 1u);
    BOOST_TEST(r.get_rows(0) == makerows(2, 42, "abc", 50, "bdef", 60, "pov"));
}

BOOST_FIXTURE_TEST_CASE(empty_row_batch, fixture)
{
    // Initial
    add_meta(r, create_meta_r1());

    // No rows, directly eof
    r.on_row_batch_start();
    auto err = r.on_row_ok_packet(create_ok_r1());
    throw_on_error(err);
    r.on_row_batch_finish();

    // Verify
    BOOST_TEST(r.is_complete());
    BOOST_TEST(r.num_resultsets() == 1u);
    BOOST_TEST(r.get_rows(0) == makerows(2));  // empty but with 2 cols
}

BOOST_FIXTURE_TEST_CASE(error_deserializing_row, fixture)
{
    add_meta(r, create_meta_r1());
    auto bad_row = create_text_row_body(42, "abc");
    bad_row.push_back(0xff);

    r.on_row_batch_start();
    auto err = r.on_row(bad_row, output_ref(), fields);
    r.on_row_batch_finish();

    BOOST_TEST(err == client_errc::extra_bytes);
}

BOOST_FIXTURE_TEST_CASE(meta_mode_minimal, fixture)
{
    exec_access(r)
        .reset(resultset_encoding::text, metadata_mode::minimal)
        .meta(create_meta_r1())
        .ok(create_ok_r1());

    BOOST_TEST(r.get_meta(0)[0].column_name() == "");
}

BOOST_FIXTURE_TEST_CASE(meta_mode_full, fixture)
{
    exec_access(r)
        .reset(resultset_encoding::text, metadata_mode::full)
        .meta(create_meta_r1())
        .ok(create_ok_r1());

    BOOST_TEST(r.get_meta(0)[0].column_name() == "ftiny");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace