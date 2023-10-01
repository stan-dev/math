//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
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

#include <boost/mysql/impl/internal/channel/channel.hpp>
#include <boost/mysql/impl/internal/network_algorithms/read_resultset_head.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/check_meta.hpp"
#include "test_common/create_diagnostics.hpp"
#include "test_unit/create_channel.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/mock_execution_processor.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::mysql::detail::channel;
using boost::mysql::detail::execution_processor;

namespace {

BOOST_AUTO_TEST_SUITE(test_read_resultset_head)

BOOST_AUTO_TEST_SUITE(detail_)  // tests the overload that can be passed an execution_processor

using netfun_maker = netfun_maker_fn<void, channel&, execution_processor&>;

struct
{
    typename netfun_maker::signature read_resultset_head;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&detail::read_resultset_head_impl),           "sync_errc"    },
    {netfun_maker::async_errinfo(&detail::async_read_resultset_head_impl), "async_errinfo"},
};

struct fixture
{
    channel chan{create_channel()};
    mock_execution_processor st;

    fixture()
    {
        // The initial request writing should have advanced this to 1 (or bigger)
        st.sequence_number() = 1;
    }

    test_stream& stream() noexcept { return get_stream(chan); }
};

BOOST_AUTO_TEST_CASE(success_meta)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()));

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_no_error();

            // We've read the response
            fix.st.num_calls().on_num_meta(1).on_meta(1).validate();
            BOOST_TEST(fix.st.is_reading_rows());
            BOOST_TEST(fix.st.sequence_number() == 3u);
            BOOST_TEST(fix.st.num_meta() == 1u);
            check_meta(fix.st.meta(), {std::make_pair(column_type::varchar, "mycol")});
        }
    }
}

BOOST_AUTO_TEST_CASE(success_several_meta_separate)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x02}))
                .add_bytes(create_coldef_frame(
                    2,
                    meta_builder().type(column_type::varchar).name("f1").build_coldef()
                ))
                .add_break()
                .add_bytes(create_coldef_frame(
                    3,
                    meta_builder().type(column_type::tinyint).name("f2").build_coldef()
                ));

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_no_error();

            // We've read the response
            fix.st.num_calls().on_num_meta(1).on_meta(2).validate();
            BOOST_TEST(fix.st.is_reading_rows());
            BOOST_TEST(fix.st.sequence_number() == 4u);
            BOOST_TEST(fix.st.num_meta() == 2u);
            check_meta(
                fix.st.meta(),
                {
                    std::make_pair(column_type::varchar, "f1"),
                    std::make_pair(column_type::tinyint, "f2"),
                }
            );
        }
    }
}

BOOST_AUTO_TEST_CASE(success_ok_packet)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(create_ok_frame(1, ok_builder().affected_rows(42).info("abc").build()));

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_no_error();

            // We've read the response
            fix.st.num_calls().on_head_ok_packet(1).validate();
            BOOST_TEST(fix.st.meta().size() == 0u);
            BOOST_TEST(fix.st.is_complete());
            BOOST_TEST(fix.st.affected_rows() == 42u);
            BOOST_TEST(fix.st.info() == "abc");
        }
    }
}

// Check that we don't attempt to read the rows even if they're available
BOOST_AUTO_TEST_CASE(success_rows_available)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_bytes(create_coldef_frame(
                    2,
                    meta_builder().type(column_type::varchar).name("f1").build_coldef()
                ))
                .add_bytes(create_text_row_message(3, "abc"));

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_no_error();

            // We've read the response but not the rows
            fix.st.num_calls().on_num_meta(1).on_meta(1).validate();
            BOOST_TEST(fix.st.is_reading_rows());
            BOOST_TEST(fix.st.sequence_number() == 3u);
        }
    }
}

// Check that we don't attempt to read the next resultset even if it's available
BOOST_AUTO_TEST_CASE(success_ok_packet_next_resultset)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_ok_frame(1, ok_builder().info("1st").more_results(true).build()))
                .add_bytes(create_ok_frame(2, ok_builder().info("2nd").build()));

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_no_error();

            // We've read the response
            fix.st.num_calls().on_head_ok_packet(1).validate();
            BOOST_TEST(fix.st.is_reading_first_subseq());
            BOOST_TEST(fix.st.info() == "1st");
        }
    }
}

// Should be a no-op
BOOST_AUTO_TEST_CASE(state_complete)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            add_ok(fix.st, ok_builder().affected_rows(42).build());

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_no_error();

            // Nothing changed
            fix.st.num_calls().on_head_ok_packet(1).validate();
            BOOST_TEST(fix.st.is_complete());
            BOOST_TEST(fix.st.affected_rows() == 42u);
        }
    }
}

// Should be a no-op
BOOST_AUTO_TEST_CASE(state_reading_rows)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            add_meta(fix.st, {meta_builder().type(column_type::bit).build_coldef()});

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_no_error();

            // Nothing changed
            fix.st.num_calls().on_num_meta(1).on_meta(1).validate();
            BOOST_TEST(fix.st.is_reading_rows());
            check_meta(fix.st.meta(), {column_type::bit});
        }
    }
}

BOOST_AUTO_TEST_CASE(error_network_error)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            // This covers: error reading the initial response and
            // error reading successive metadata packets
            for (std::size_t i = 0; i <= 2; ++i)
            {
                BOOST_TEST_CONTEXT(i)
                {
                    fixture fix;
                    fix.stream()
                        .add_bytes(create_frame(1, {0x02}))
                        .add_break()
                        .add_bytes(create_coldef_frame(
                            2,
                            meta_builder().type(column_type::varchar).name("f1").build_coldef()
                        ))
                        .add_break()
                        .add_bytes(create_coldef_frame(
                            3,
                            meta_builder().type(column_type::tinyint).name("f2").build_coldef()
                        ))
                        .set_fail_count(fail_count(i, client_errc::server_unsupported));

                    // Call the function
                    fns.read_resultset_head(fix.chan, fix.st)
                        .validate_error_exact(client_errc::server_unsupported);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(error_metadata_packets_seqnum_mismatch)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x02}))
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()))
                .add_bytes(create_coldef_frame(4, meta_builder().type(column_type::tinyint).build_coldef()));

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st)
                .validate_error_exact(client_errc::sequence_number_mismatch);
        }
    }
}

// All cases where the deserialization of the execution_response
// yields an error are handled uniformly, so it's enough with this test
BOOST_AUTO_TEST_CASE(error_deserialize_execution_response)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(err_builder()
                                       .seqnum(1)
                                       .code(common_server_errc::er_bad_db_error)
                                       .message("no_db")
                                       .build_frame());

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st)
                .validate_error_exact(common_server_errc::er_bad_db_error, "no_db");
        }
    }
}

BOOST_AUTO_TEST_CASE(error_deserialize_metadata)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_bytes(create_frame(2, {0x08, 0x03}));  // bad coldef

            // Call the function
            fns.read_resultset_head(fix.chan, fix.st).validate_error_exact(client_errc::incomplete_message);
        }
    }
}

// The execution processor signals an error on head packet (e.g. meta mismatch)
BOOST_AUTO_TEST_CASE(error_on_head_ok_packet)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.st.set_fail_count(
                fail_count(0, client_errc::metadata_check_failed),
                create_client_diag("some message")
            );

            fix.stream().add_bytes(create_ok_frame(1, ok_builder().affected_rows(42).info("abc").build()));

            fns.read_resultset_head(fix.chan, fix.st)
                .validate_error_exact_client(client_errc::metadata_check_failed, "some message");
            fix.st.num_calls().on_head_ok_packet(1).validate();
        }
    }
}

BOOST_AUTO_TEST_CASE(error_on_meta)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.st.set_fail_count(
                fail_count(0, client_errc::metadata_check_failed),
                create_client_diag("some message")
            );

            fix.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()));

            fns.read_resultset_head(fix.chan, fix.st)
                .validate_error_exact_client(client_errc::metadata_check_failed, "some message");
            fix.st.num_calls().on_num_meta(1).on_meta(1).validate();
        }
    }
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace