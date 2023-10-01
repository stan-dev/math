//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/metadata_mode.hpp>

#include <boost/mysql/detail/any_execution_request.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>
#include <boost/mysql/impl/internal/network_algorithms/start_execution.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/check_meta.hpp"
#include "test_common/create_basic.hpp"
#include "test_unit/create_channel.hpp"
#include "test_unit/create_coldef_frame.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_statement.hpp"
#include "test_unit/mock_execution_processor.hpp"
#include "test_unit/printing.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::mysql::detail::any_execution_request;
using boost::mysql::detail::channel;
using boost::mysql::detail::execution_processor;
using boost::mysql::detail::resultset_encoding;

BOOST_AUTO_TEST_SUITE(test_start_execution)

using netfun_maker = netfun_maker_fn<void, channel&, const any_execution_request&, execution_processor&>;

struct
{
    typename netfun_maker::signature start_execution;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&boost::mysql::detail::start_execution_impl),           "sync" },
    {netfun_maker::async_errinfo(&boost::mysql::detail::async_start_execution_impl), "async"}
};

struct fixture
{
    channel chan{create_channel()};
    mock_execution_processor st;

    fixture() { chan.shared_sequence_number() = 42; }

    test_stream& stream() { return get_stream(chan); }
};

BOOST_AUTO_TEST_CASE(text_query)
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
            fns.start_execution(fix.chan, any_execution_request("SELECT 1"), fix.st).validate_no_error();

            // We've written the request message
            auto expected_msg = create_frame(0, {0x03, 0x53, 0x45, 0x4c, 0x45, 0x43, 0x54, 0x20, 0x31});
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), expected_msg);
            BOOST_TEST(fix.chan.shared_sequence_number() == 42u);  // unused

            // We've read the response
            BOOST_TEST(fix.st.encoding() == resultset_encoding::text);
            BOOST_TEST(fix.st.sequence_number() == 3u);
            BOOST_TEST(fix.st.is_reading_rows());
            check_meta(fix.st.meta(), {column_type::varchar});

            // Validate mock calls
            fix.st.num_calls().reset(1).on_num_meta(1).on_meta(1).validate();
        }
    }
}

BOOST_AUTO_TEST_CASE(prepared_statement)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()));
            auto stmt = statement_builder().id(1).num_params(2).build();
            const auto params = make_fv_arr("test", nullptr);

            // Call the function
            fns.start_execution(fix.chan, any_execution_request(stmt, params), fix.st).validate_no_error();

            // We've written the request message
            constexpr std::uint8_t body[] = {
                0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02,
                0x01, 0xfe, 0x00, 0x06, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74,
            };
            auto expected_msg = create_frame(0, body);
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), expected_msg);
            BOOST_TEST(fix.chan.shared_sequence_number() == 42u);  // unused

            // We've read the response
            BOOST_TEST(fix.st.encoding() == resultset_encoding::binary);
            BOOST_TEST(fix.st.sequence_number() == 3u);
            BOOST_TEST(fix.st.is_reading_rows());
            check_meta(fix.st.meta(), {column_type::varchar});

            // Validate mock calls
            fix.st.num_calls().reset(1).on_num_meta(1).on_meta(1).validate();
        }
    }
}

BOOST_AUTO_TEST_CASE(error_num_params)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream()
                .add_bytes(create_frame(1, {0x01}))
                .add_bytes(create_coldef_frame(2, meta_builder().type(column_type::varchar).build_coldef()));
            auto stmt = statement_builder().id(1).num_params(2).build();
            const auto params = make_fv_arr("test", nullptr, 42);  // too many params

            // Call the function
            fns.start_execution(fix.chan, any_execution_request(stmt, params), fix.st)
                .validate_error_exact(client_errc::wrong_num_params);

            // We didn't write any message and didn't modify the processor
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), std::vector<std::uint8_t>());
            fix.st.num_calls().validate();
        }
    }
}

// This covers errors in both writing the request and calling read_resultset_head
BOOST_AUTO_TEST_CASE(error_network_error)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            for (std::size_t i = 0; i <= 1; ++i)
            {
                BOOST_TEST_CONTEXT(i)
                {
                    fixture fix;
                    fix.stream()
                        .add_bytes(create_frame(1, {0x01}))
                        .add_bytes(
                            create_coldef_frame(2, meta_builder().type(column_type::bit).build_coldef())
                        )
                        .set_fail_count(fail_count(i, client_errc::server_unsupported));

                    // Call the function
                    fns.start_execution(fix.chan, any_execution_request("SELECT 1"), fix.st)
                        .validate_error_exact(client_errc::server_unsupported);

                    // Num calls validation
                    fix.st.num_calls().reset(1).validate();
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
