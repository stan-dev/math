//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>

#include <boost/mysql/impl/internal/sansio/connection_state_data.hpp>
#include <boost/mysql/impl/internal/sansio/read_some_rows.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>

#include "test_common/buffer_concat.hpp"
#include "test_common/create_diagnostics.hpp"
#include "test_unit/algo_test.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/mock_execution_processor.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::span;
using boost::mysql::detail::output_ref;

BOOST_AUTO_TEST_SUITE(test_read_some_rows)

struct fixture : algo_fixture_base
{
    using row1 = std::tuple<int>;

    mock_execution_processor proc;
    std::array<row1, 3> storage;
    detail::read_some_rows_algo algo{
        diag,
        {&proc, ref()}
    };

    output_ref ref() noexcept { return output_ref(span<row1>(storage), 0); }

    fixture()
    {
        // Prepare the processor, such that it's ready to read rows
        add_meta(
            proc,
            {meta_builder().type(column_type::varchar).name("fvarchar").nullable(false).build_coldef()}
        );
        proc.sequence_number() = 42;
        st.backslash_escapes = false;
    }

    void validate_refs(std::size_t num_rows)
    {
        BOOST_TEST_REQUIRE(proc.refs().size() == num_rows);
        for (std::size_t i = 0; i < num_rows; ++i)
            BOOST_TEST(proc.refs()[i].offset() == i);
    }

    std::size_t result() const { return algo.result(st); }
};

BOOST_AUTO_TEST_CASE(eof)
{
    // Setup
    fixture fix;

    // Run the test
    algo_test()
        .expect_read(
            create_eof_frame(42, ok_builder().affected_rows(1).info("1st").more_results(true).build())
        )
        .check(fix);

    BOOST_TEST(fix.result() == 0u);  // num read rows
    BOOST_TEST(fix.proc.is_reading_head());
    BOOST_TEST(fix.proc.affected_rows() == 1u);
    BOOST_TEST(fix.proc.info() == "1st");
    BOOST_TEST(fix.st.backslash_escapes);
}

BOOST_AUTO_TEST_CASE(eof_no_backslash_escapes)
{
    // Setup
    fixture fix;

    // Run the test
    algo_test()
        .expect_read(create_eof_frame(42, ok_builder().no_backslash_escapes(true).more_results(true).build()))
        .check(fix);

    BOOST_TEST(fix.result() == 0u);  // num read rows
    BOOST_TEST(fix.proc.is_reading_head());
    BOOST_TEST(!fix.st.backslash_escapes);
}

BOOST_AUTO_TEST_CASE(batch_with_rows)
{
    // Setup
    fixture fix;

    // Run the algo. Single, long read that yields two rows
    algo_test()
        .expect_read(buffer_builder()
                         .add(create_text_row_message(42, "abc"))
                         .add(create_text_row_message(43, "von"))
                         .build())
        .check(fix);

    // Validate
    BOOST_TEST(fix.result() == 2u);  // num read rows
    BOOST_TEST(fix.proc.is_reading_rows());
    fix.validate_refs(2);
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(1)
        .on_row(2)
        .on_row_batch_finish(1)
        .validate();
}

BOOST_AUTO_TEST_CASE(batch_with_rows_eof)
{
    // Setup
    fixture fix;

    // Run the algo. Single, long read that yields rows and eof
    algo_test()
        .expect_read(
            buffer_builder()
                .add(create_text_row_message(42, "abc"))
                .add(create_text_row_message(43, "von"))
                .add(
                    create_eof_frame(44, ok_builder().affected_rows(1).info("1st").more_results(true).build())
                )
                .build()
        )
        .check(fix);

    // Validate
    BOOST_TEST(fix.result() == 2u);  // num read rows
    BOOST_TEST_REQUIRE(fix.proc.is_reading_head());
    BOOST_TEST(fix.proc.affected_rows() == 1u);
    BOOST_TEST(fix.proc.info() == "1st");
    fix.validate_refs(2);
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(1)
        .on_row(2)
        .on_row_ok_packet(1)
        .on_row_batch_finish(1)
        .validate();
}

// Regression check: don't attempt to continue reading after the 1st EOF for multi-result
BOOST_AUTO_TEST_CASE(batch_with_rows_eof_multiresult)
{
    // Setup
    fixture fix;

    // Run the algo. Single, long read that yields the next resultset's OK packet
    algo_test()
        .expect_read(
            buffer_builder()
                .add(create_text_row_message(42, "abc"))
                .add(
                    create_eof_frame(43, ok_builder().affected_rows(1).info("1st").more_results(true).build())
                )
                .add(create_ok_frame(44, ok_builder().info("2nd").build()))
                .build()
        )
        .check(fix);

    // Validate
    BOOST_TEST(fix.result() == 1u);  // num read rows
    BOOST_TEST_REQUIRE(fix.proc.is_reading_head());
    BOOST_TEST(fix.proc.affected_rows() == 1u);
    BOOST_TEST(fix.proc.info() == "1st");
    fix.validate_refs(1);
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(1)
        .on_row(1)
        .on_row_ok_packet(1)
        .on_row_batch_finish(1)
        .validate();
}

BOOST_AUTO_TEST_CASE(batch_with_rows_out_of_span_space)
{
    // Setup
    fixture fix;

    // Run the algo. Single, long read that yields 4 rows.
    // We have only space for 3
    algo_test()
        .expect_read(buffer_builder()
                         .add(create_text_row_message(42, "aaa"))
                         .add(create_text_row_message(43, "bbb"))
                         .add(create_text_row_message(44, "ccc"))
                         .add(create_text_row_message(45, "ddd"))
                         .build())
        .check(fix);

    // Validate
    BOOST_TEST(fix.result() == 3u);  // num read rows
    fix.validate_refs(3);
    BOOST_TEST(fix.proc.is_reading_rows());
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(1)
        .on_row(3)
        .on_row_batch_finish(1)
        .validate();
}

BOOST_AUTO_TEST_CASE(successive_calls_keep_parsing_state)
{
    // Setup
    fixture fix;

    // Run the algo
    auto eof = create_eof_frame(44, ok_builder().affected_rows(1).info("1st").build());
    algo_test()
        .expect_read(buffer_builder()
                         .add(create_text_row_message(42, "aaa"))
                         .add(create_text_row_message(43, "bbb"))
                         .add(boost::span<const std::uint8_t>(eof).subspan(0, 6))  // OK partially received
                         .build())
        .check(fix);

    // Validate
    BOOST_TEST(fix.result() == 2u);  // num read rows
    fix.validate_refs(2);
    BOOST_TEST(fix.proc.is_reading_rows());
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(1)
        .on_row(2)
        .on_row_batch_finish(1)
        .validate();

    // Run the algo again
    fix.algo.reset();
    fix.diag = create_server_diag("Diagnostics not cleared");
    algo_test()
        .expect_read(buffer_builder().add(boost::span<const std::uint8_t>(eof).subspan(6)).build())
        .check(fix);

    // Validate
    BOOST_TEST(fix.result() == 0u);  // num read rows
    BOOST_TEST(fix.proc.is_complete());
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(2)
        .on_row(2)
        .on_row_batch_finish(2)
        .on_row_ok_packet(1)
        .validate();
    BOOST_TEST(fix.proc.affected_rows() == 1u);
    BOOST_TEST(fix.proc.info() == "1st");
}

// read_some_rows is a no-op if !st.should_read_rows()
BOOST_AUTO_TEST_CASE(state_complete)
{
    // Setup
    fixture fix;
    add_ok(fix.proc, ok_builder().affected_rows(20).build());

    // Run the algo
    algo_test().check(fix);

    // Validate
    BOOST_TEST(fix.result() == 0u);  // num read rows
    BOOST_TEST(fix.proc.is_complete());
    fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row_ok_packet(1).validate();
}

BOOST_AUTO_TEST_CASE(state_reading_head)
{
    // Setup
    fixture fix;
    add_ok(fix.proc, ok_builder().affected_rows(42).more_results(true).build());

    // Run the algo
    algo_test().check(fix);

    // Validate
    BOOST_TEST(fix.result() == 0u);  // num read rows
    BOOST_TEST(fix.proc.is_reading_head());
    fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row_ok_packet(1).validate();
}

BOOST_AUTO_TEST_CASE(error_network_error)
{
    algo_test().expect_read(create_text_row_message(42, "aaa")).check_network_errors<fixture>();
}

BOOST_AUTO_TEST_CASE(error_seqnum_mismatch_successive_messages)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(buffer_builder()
                         .add(create_text_row_message(42, "abc"))
                         .add(create_text_row_message(45, "von"))  // seqnum mismatch here
                         .build())
        .check(fix, client_errc::sequence_number_mismatch);
}

BOOST_AUTO_TEST_CASE(error_on_row)
{
    // Setup
    fixture fix;

    // Mock a failure
    fix.proc.set_fail_count(fail_count(0, client_errc::static_row_parsing_error));

    // Run the algo
    algo_test()
        .expect_read(create_text_row_message(42, 10))
        .check(fix, client_errc::static_row_parsing_error);

    // Validate
    fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row(1).on_row_batch_start(1).validate();
}

BOOST_AUTO_TEST_CASE(error_on_row_ok_packet)
{
    // Setup
    fixture fix;

    // Mock a failure
    fix.proc.set_fail_count(fail_count(0, client_errc::num_resultsets_mismatch));

    // Run the algo
    algo_test()
        .expect_read(create_eof_frame(42, ok_builder().build()))
        .check(fix, client_errc::num_resultsets_mismatch);

    // Validate
    fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row_ok_packet(1).on_row_batch_start(1).validate();
}

// deserialize_row_message covers cases like getting an error packet, deserialization errors, etc.
BOOST_AUTO_TEST_CASE(error_deserialize_row_message)
{
    // Setup
    fixture fix;

    // Run the algo
    algo_test()
        .expect_read(
            err_builder().seqnum(42).code(common_server_errc::er_alter_info).message("abc").build_frame()
        )
        .check(fix, common_server_errc::er_alter_info, create_server_diag("abc"));

    // Validate
    fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row_batch_start(1).validate();
}

BOOST_AUTO_TEST_CASE(reset)
{
    // Setup
    fixture fix;

    // Run the algo. Read a row
    algo_test().expect_read(create_text_row_message(42, "abc")).check(fix);

    // Validate
    BOOST_TEST(fix.result() == 1u);  // num read rows
    BOOST_TEST(fix.proc.is_reading_rows());
    fix.validate_refs(1);
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(1)
        .on_row(1)
        .on_row_batch_finish(1)
        .validate();

    // Reset
    fix.algo.reset();

    // Run the algo again. Read an OK packet
    algo_test().expect_read(create_eof_frame(43, ok_builder().build())).check(fix);

    // Check
    BOOST_TEST(fix.result() == 0u);  // num read rows
    BOOST_TEST(fix.proc.is_complete());
    fix.proc.num_calls()
        .on_num_meta(1)
        .on_meta(1)
        .on_row_batch_start(2)
        .on_row(1)
        .on_row_batch_finish(2)
        .on_row_ok_packet(1)
        .validate();
}

BOOST_AUTO_TEST_SUITE_END()
