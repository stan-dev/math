//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/sansio/connection_state_data.hpp>
#include <boost/mysql/impl/internal/sansio/read_resultset_head.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/buffer_concat.hpp"
#include "test_common/check_meta.hpp"
#include "test_unit/algo_test.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/mock_execution_processor.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_read_resultset_head)

struct fixture : algo_fixture_base
{
    mock_execution_processor proc;
    detail::read_resultset_head_algo algo{diag, {&proc}};

    fixture()
    {
        // The initial request writing should have advanced this to 1 (or bigger)
        proc.sequence_number() = 1;
    }
};

BOOST_AUTO_TEST_CASE(success_meta)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(create_frame(1, {0x01}))  // 1 metadata follows
        .expect_read(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .check(fix);

    // Verify
    fix.proc.num_calls().on_num_meta(1).on_meta(1).validate();
    BOOST_TEST(fix.proc.is_reading_rows());
    BOOST_TEST(fix.proc.sequence_number() == 3u);
    BOOST_TEST(fix.proc.num_meta() == 1u);
    check_meta(fix.proc.meta(), {std::make_pair(column_type::varchar, "mycol")});
    BOOST_TEST(fix.st.backslash_escapes);
}

BOOST_AUTO_TEST_CASE(success_ok_packet)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(create_ok_frame(1, ok_builder().affected_rows(42).info("abc").build()))
        .check(fix);

    // Verify
    fix.proc.num_calls().on_head_ok_packet(1).validate();
    BOOST_TEST(fix.proc.meta().size() == 0u);
    BOOST_TEST(fix.proc.is_complete());
    BOOST_TEST(fix.proc.affected_rows() == 42u);
    BOOST_TEST(fix.proc.info() == "abc");
    BOOST_TEST(fix.st.backslash_escapes);
}

BOOST_AUTO_TEST_CASE(success_ok_packet_no_backslash_escapes)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test().expect_read(create_ok_frame(1, ok_builder().no_backslash_escapes(true).build())).check(fix);

    // Verify
    fix.proc.num_calls().on_head_ok_packet(1).validate();
    BOOST_TEST(!fix.st.backslash_escapes);
}

// Check that we don't attempt to read the rows even if they're available
BOOST_AUTO_TEST_CASE(success_rows_available)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(create_frame(1, {0x01}))  // 1 metadata follows
        .expect_read(buffer_builder()
                         .add(create_coldef_frame(
                             2,
                             meta_builder().type(column_type::varchar).name("f1").build_coldef()
                         ))
                         .add(create_text_row_message(3, "abc"))
                         .build())
        .check(fix);

    // We've read the response but not the rows
    fix.proc.num_calls().on_num_meta(1).on_meta(1).validate();
    BOOST_TEST(fix.proc.is_reading_rows());
    BOOST_TEST(fix.proc.sequence_number() == 3u);
}

// Check that we don't attempt to read the next resultset even if it's available
BOOST_AUTO_TEST_CASE(success_ok_packet_next_resultset)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(buffer_builder()
                         .add(create_ok_frame(1, ok_builder().info("1st").more_results(true).build()))
                         .add(create_ok_frame(2, ok_builder().info("2nd").build()))
                         .build())
        .check(fix);

    // Verify
    fix.proc.num_calls().on_head_ok_packet(1).validate();
    BOOST_TEST(fix.proc.is_reading_first_subseq());
    BOOST_TEST(fix.proc.info() == "1st");
}

BOOST_AUTO_TEST_CASE(state_complete)
{
    // Setup
    fixture fix;
    add_ok(fix.proc, ok_builder().affected_rows(42).build());

    // Should be a no-op
    algo_test().check(fix);

    // Nothing changed
    fix.proc.num_calls().on_head_ok_packet(1).validate();
    BOOST_TEST(fix.proc.is_complete());
    BOOST_TEST(fix.proc.affected_rows() == 42u);
}

BOOST_AUTO_TEST_CASE(state_reading_rows)
{
    // Setup
    fixture fix;
    add_meta(fix.proc, {meta_builder().type(column_type::bit).build_coldef()});

    // Should be a no-op
    algo_test().check(fix);

    // Nothing changed
    fix.proc.num_calls().on_num_meta(1).on_meta(1).validate();
    BOOST_TEST(fix.proc.is_reading_rows());
    check_meta(fix.proc.meta(), {column_type::bit});
}

BOOST_AUTO_TEST_CASE(error_network_error)
{
    // This covers testing for network errors for all the reads we perform
    algo_test()
        .expect_read(create_frame(1, {0x02}))
        .expect_read(
            create_coldef_frame(2, meta_builder().type(column_type::varchar).name("f1").build_coldef())
        )
        .expect_read(
            create_coldef_frame(3, meta_builder().type(column_type::tinyint).name("f2").build_coldef())
        )
        .check_network_errors<fixture>();
}

// All cases where the deserialization of the execution_response
// yields an error are handled uniformly, so it's enough with this test
BOOST_AUTO_TEST_CASE(error_deserialize_execution_response)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(
            err_builder().seqnum(1).code(common_server_errc::er_bad_db_error).message("no_db").build_frame()
        )
        .check(fix, common_server_errc::er_bad_db_error, create_server_diag("no_db"));
}

BOOST_AUTO_TEST_CASE(error_deserialize_metadata)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(create_frame(1, {0x01}))
        .expect_read(create_frame(2, {0x08, 0x03}))  // bad coldef
        .check(fix, client_errc::incomplete_message);
}

// The execution processor signals an error on head packet (e.g. meta mismatch)
BOOST_AUTO_TEST_CASE(error_on_head_ok_packet)
{
    // Setup
    fixture fix;
    fix.proc.set_fail_count(
        fail_count(0, client_errc::metadata_check_failed),
        create_client_diag("some message")
    );

    // Run the algo
    algo_test()
        .expect_read(create_ok_frame(1, ok_builder().affected_rows(42).info("abc").build()))
        .check(fix, client_errc::metadata_check_failed, create_client_diag("some message"));

    // Verify
    fix.proc.num_calls().on_head_ok_packet(1).validate();
}

BOOST_AUTO_TEST_CASE(error_on_meta)
{
    // Setup
    fixture fix;
    fix.proc.set_fail_count(
        fail_count(0, client_errc::metadata_check_failed),
        create_client_diag("some message")
    );

    // Run the algo
    algo_test()
        .expect_read(create_frame(1, {0x01}))
        .expect_read(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .check(fix, client_errc::metadata_check_failed, create_client_diag("some message"));

    // Verify
    fix.proc.num_calls().on_num_meta(1).on_meta(1).validate();
}

BOOST_AUTO_TEST_CASE(reset)
{
    // Setup
    fixture fix;

    // Run the algo once
    algo_test()
        .expect_read(create_frame(1, {0x01}))  // 1 metadata follows
        .expect_read(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
        .check(fix);
    fix.proc.num_calls().on_num_meta(1).on_meta(1).validate();

    // Reset. Place the processor into a state where we can read head again
    fix.algo.reset();
    add_ok(fix.proc, ok_builder().more_results(true).build());

    // Run it again
    algo_test().expect_read(create_ok_frame(3, ok_builder().build())).check(fix);
    fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row_ok_packet(1).on_head_ok_packet(1).validate();
}

BOOST_AUTO_TEST_SUITE_END()
