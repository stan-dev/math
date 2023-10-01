//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>
#include <boost/mysql/impl/internal/network_algorithms/read_some_rows.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>

#include "test_unit/create_channel.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/mock_execution_processor.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::span;
using boost::mysql::detail::channel;
using boost::mysql::detail::execution_processor;
using boost::mysql::detail::output_ref;

BOOST_AUTO_TEST_SUITE(test_read_some_rows)

using row1 = std::tuple<int>;

using netfun_maker = netfun_maker_fn<std::size_t, channel&, execution_processor&, const output_ref&>;

struct
{
    typename netfun_maker::signature read_some_rows_impl;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&detail::read_some_rows_impl),           "sync" },
    {netfun_maker::async_errinfo(&detail::async_read_some_rows_impl), "async"},
};

struct fixture
{
    mock_execution_processor proc;
    channel chan{create_channel()};
    std::array<row1, 3> storage;

    fixture()
    {
        // Prepare the processor, such that it's ready to read rows
        add_meta(
            proc,
            {meta_builder().type(column_type::varchar).name("fvarchar").nullable(false).build_coldef()}
        );
        proc.sequence_number() = 42;
    }

    void validate_refs(std::size_t num_rows)
    {
        BOOST_TEST_REQUIRE(proc.refs().size() == num_rows);
        for (std::size_t i = 0; i < num_rows; ++i)
            BOOST_TEST(proc.refs()[i].offset() == i);
    }

    test_stream& stream() noexcept { return get_stream(chan); }
    output_ref ref() noexcept { return output_ref(span<row1>(storage), 0); }
};

BOOST_AUTO_TEST_CASE(eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(
                create_eof_frame(42, ok_builder().affected_rows(1).info("1st").more_results(true).build())
            );

            std::size_t num_rows = fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref()).get();
            BOOST_TEST(num_rows == 0u);
            BOOST_TEST(fix.proc.is_reading_head());
            BOOST_TEST(fix.proc.affected_rows() == 1u);
            BOOST_TEST(fix.proc.info() == "1st");
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
            fix.proc.num_calls()
                .on_num_meta(1)
                .on_meta(1)
                .on_row_batch_start(1)
                .on_row_ok_packet(1)
                .on_row_batch_finish(1)
                .validate();
        }
    }
}

BOOST_AUTO_TEST_CASE(batch_with_rows)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_text_row_message(42, "abc"))
                .add_bytes(create_text_row_message(43, "von"))
                .add_break()
                .add_bytes(create_text_row_message(44, "other"));  // only a single read should be issued

            std::size_t num_rows = fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref()).get();
            BOOST_TEST(num_rows == 2u);
            BOOST_TEST(fix.proc.is_reading_rows());
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
            fix.validate_refs(2);
            fix.proc.num_calls()
                .on_num_meta(1)
                .on_meta(1)
                .on_row_batch_start(1)
                .on_row(2)
                .on_row_batch_finish(1)
                .validate();
        }
    }
}

BOOST_AUTO_TEST_CASE(batch_with_rows_eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_text_row_message(42, "abc"))
                .add_bytes(create_text_row_message(43, "von"))
                .add_bytes(
                    create_eof_frame(44, ok_builder().affected_rows(1).info("1st").more_results(true).build())
                );

            std::size_t num_rows = fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref()).get();
            BOOST_TEST(num_rows == 2u);
            BOOST_TEST_REQUIRE(fix.proc.is_reading_head());
            BOOST_TEST(fix.proc.affected_rows() == 1u);
            BOOST_TEST(fix.proc.info() == "1st");
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
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
    }
}

// Regression check: don't attempt to continue reading after the 1st EOF for multi-result
BOOST_AUTO_TEST_CASE(batch_with_rows_eof_multiresult)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_text_row_message(42, "abc"))
                .add_bytes(
                    create_eof_frame(43, ok_builder().affected_rows(1).info("1st").more_results(true).build())
                )
                .add_bytes(create_ok_frame(44, ok_builder().info("2nd").build()));

            std::size_t num_rows = fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref()).get();
            BOOST_TEST(num_rows == 1u);
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
    }
}

BOOST_AUTO_TEST_CASE(batch_with_rows_out_of_span_space)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_text_row_message(42, "aaa"))
                .add_bytes(create_text_row_message(43, "bbb"))
                .add_bytes(create_text_row_message(44, "ccc"))
                .add_bytes(create_text_row_message(45, "ddd"));

            // We only have space for 3
            std::size_t num_rows = fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref()).get();
            BOOST_TEST(num_rows == 3u);
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
    }
}

// read_some_rows is a no-op if !st.should_read_rows()
BOOST_AUTO_TEST_CASE(state_complete)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            add_ok(fix.proc, ok_builder().affected_rows(20).build());

            std::size_t num_rows = fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref()).get();
            BOOST_TEST(num_rows == 0u);
            BOOST_TEST(fix.proc.is_complete());
            fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row_ok_packet(1).validate();
        }
    }
}

BOOST_AUTO_TEST_CASE(state_reading_head)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            add_ok(fix.proc, ok_builder().affected_rows(42).more_results(true).build());

            std::size_t num_rows = fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref()).get();
            BOOST_TEST(num_rows == 0u);
            BOOST_TEST(fix.proc.is_reading_head());
            fix.proc.num_calls().on_num_meta(1).on_meta(1).on_row_ok_packet(1).validate();
        }
    }
}

BOOST_AUTO_TEST_CASE(error_network_error)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            for (std::size_t i = 0; i <= 1; ++i)
            {
                BOOST_TEST_CONTEXT("i=" << i)
                {
                    fixture fix;
                    fix.stream()
                        .add_bytes(create_text_row_message(42, "abc"))
                        .add_break()
                        .add_bytes(create_eof_frame(43, ok_builder().affected_rows(1).info("1st").build()))
                        .set_fail_count(fail_count(i, client_errc::wrong_num_params));

                    fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref())
                        .validate_error_exact(client_errc::wrong_num_params);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(error_on_row)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(create_text_row_message(42, 10));

            // Mock a failure
            fix.proc.set_fail_count(fail_count(0, client_errc::static_row_parsing_error));

            // Call the function
            fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref())
                .validate_error_exact(client_errc::static_row_parsing_error);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_on_row_ok_packet)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(create_eof_frame(42, ok_builder().build()));

            // Mock a failure
            fix.proc.set_fail_count(fail_count(0, client_errc::num_resultsets_mismatch));

            // Call the function
            fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref())
                .validate_error_exact(client_errc::num_resultsets_mismatch);
        }
    }
}

// deserialize_row_message covers cases like getting an error packet, deserialization errors, etc.
BOOST_AUTO_TEST_CASE(error_deserialize_row_message)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(
                err_builder().seqnum(42).code(common_server_errc::er_alter_info).build_frame()
            );

            // Call the function
            fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref())
                .validate_error_exact(common_server_errc::er_alter_info);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_seqnum_mismatch)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(create_eof_frame(50, ok_builder().build()));

            // Call the function
            fns.read_some_rows_impl(fix.chan, fix.proc, fix.ref())
                .validate_error_exact(client_errc::sequence_number_mismatch);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
