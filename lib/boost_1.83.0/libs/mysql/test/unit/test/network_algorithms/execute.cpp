//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>

#include <boost/mysql/detail/any_execution_request.hpp>
#include <boost/mysql/detail/execution_processor/execution_processor.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>
#include <boost/mysql/impl/internal/network_algorithms/execute.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/buffer_concat.hpp"
#include "test_common/check_meta.hpp"
#include "test_unit/create_channel.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/mock_execution_processor.hpp"
#include "test_unit/printing.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::mysql::detail::any_execution_request;
using boost::mysql::detail::channel;
using boost::mysql::detail::execution_processor;
using boost::mysql::detail::resultset_encoding;

BOOST_AUTO_TEST_SUITE(test_execute)

using netfun_maker = netfun_maker_fn<void, channel&, const any_execution_request&, execution_processor&>;

struct
{
    typename netfun_maker::signature execute;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&detail::execute_impl),           "sync" },
    {netfun_maker::async_errinfo(&detail::async_execute_impl), "async"}
};

struct fixture
{
    mock_execution_processor proc;
    channel chan{create_channel()};

    test_stream& stream() noexcept { return get_stream(chan); }
};

// The serialized form of a SELECT 1 query request
constexpr std::uint8_t serialized_select_1[] = {0x03, 0x53, 0x45, 0x4c, 0x45, 0x43, 0x54, 0x20, 0x31};

BOOST_AUTO_TEST_CASE(eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(create_ok_frame(1, ok_builder().affected_rows(60u).info("abc").build()));

            // Call the function
            fns.execute(fix.chan, any_execution_request("SELECT 1"), fix.proc).validate_no_error();

            // We've written the exeuction request
            auto expected_msg = create_frame(0, serialized_select_1);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), expected_msg);

            // We've read into the processor
            fix.proc.num_calls().reset(1).on_head_ok_packet(1).validate();
            BOOST_TEST(fix.proc.encoding() == resultset_encoding::text);
            BOOST_TEST(fix.proc.affected_rows() == 60u);
            BOOST_TEST(fix.proc.info() == "abc");
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
        }
    }
}

BOOST_AUTO_TEST_CASE(single_batch)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x01}))  // OK, 1 column
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::bigint).build_coldef()))
                .add_bytes(create_text_row_message(3, 42))  // row 1
                .add_bytes(create_text_row_message(4, 43))  // row 2
                .add_bytes(create_eof_frame(5, ok_builder().affected_rows(10u).info("1st").build()))
                .add_bytes(create_ok_frame(1, ok_builder().info("2nd").build()));  // don't read any further

            // Call the function
            fns.execute(fix.chan, any_execution_request("SELECT 1"), fix.proc).validate_no_error();

            // We've written the exeuction request
            auto expected_msg = create_frame(0, serialized_select_1);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), expected_msg);

            // We've read the results
            fix.proc.num_calls()
                .reset(1)
                .on_num_meta(1)
                .on_meta(1)
                .on_row_batch_start(1)
                .on_row(2)
                .on_row_batch_finish(1)
                .on_row_ok_packet(1)
                .validate();
            BOOST_TEST(fix.proc.encoding() == resultset_encoding::text);
            BOOST_TEST(fix.proc.num_meta() == 1u);
            check_meta(fix.proc.meta(), {column_type::bigint});
            BOOST_TEST(fix.proc.affected_rows() == 10u);
            BOOST_TEST(fix.proc.info() == "1st");
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
        }
    }
}

BOOST_AUTO_TEST_CASE(multiple_batches)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x01}))  // OK, 1 column
                .add_break()
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::tinyint).build_coldef()))
                .add_break()
                .add_bytes(create_text_row_message(3, 42))  // row 1
                .add_break()
                .add_bytes(create_text_row_message(4, 43))  // row 2
                .add_break()
                .add_bytes(create_eof_frame(5, ok_builder().affected_rows(10u).info("1st").build()));

            // Call the function
            fns.execute(fix.chan, any_execution_request("SELECT 1"), fix.proc).validate_no_error();

            // We've written the exeuction request
            auto expected_msg = create_frame(0, serialized_select_1);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), expected_msg);

            // We've read the results
            fix.proc.num_calls()
                .reset(1)
                .on_num_meta(1)
                .on_meta(1)
                .on_row_batch_start(3)
                .on_row(2)
                .on_row_batch_finish(3)
                .on_row_ok_packet(1)
                .validate();
            BOOST_TEST(fix.proc.encoding() == resultset_encoding::text);
            BOOST_TEST(fix.proc.num_meta() == 1u);
            check_meta(fix.proc.meta(), {column_type::tinyint});
            BOOST_TEST(fix.proc.affected_rows() == 10u);
            BOOST_TEST(fix.proc.info() == "1st");
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
        }
    }
}

// Tests error on write, while reading head and while reading rows (error spotcheck)
BOOST_AUTO_TEST_CASE(error_network_error)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            for (std::size_t i = 0; i <= 2; ++i)
            {
                BOOST_TEST_CONTEXT("i=" << i)
                {
                    fixture fix;
                    fix.stream()
                        .add_bytes(create_frame(1, {0x01}))  // OK, 1 column
                        .add_bytes(
                            create_coldef_frame(2, meta_builder().type(column_type::tinyint).build_coldef())
                        )
                        .add_break()
                        .add_bytes(create_text_row_message(3, 42))
                        .add_bytes(create_eof_frame(4, ok_builder().info("1st").build()))
                        .set_fail_count(fail_count(i, client_errc::wrong_num_params));

                    // Call the function
                    fns.execute(fix.chan, any_execution_request("SELECT 1"), fix.proc)
                        .validate_error_exact(client_errc::wrong_num_params);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
