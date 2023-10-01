//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>

#include <boost/mysql/detail/execution_processor/execution_state_impl.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>
#include <boost/mysql/impl/internal/network_algorithms/read_some_rows_dynamic.hpp>

#include <boost/test/unit_test.hpp>

#include "test_unit/create_channel.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::mysql::detail::channel;
using boost::mysql::detail::execution_state_impl;

BOOST_AUTO_TEST_SUITE(test_read_some_rows_dynamic)

using netfun_maker = netfun_maker_fn<rows_view, channel&, execution_state_impl&>;

struct
{
    typename netfun_maker::signature read_some_rows_dynamic;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&detail::read_some_rows_dynamic_impl),           "sync" },
    {netfun_maker::async_errinfo(&detail::async_read_some_rows_dynamic_impl), "async"},
};

struct fixture
{
    execution_state_impl st;
    channel chan{create_channel()};

    fixture()
    {
        // Prepare the state, such that it's ready to read rows
        add_meta(st, {meta_builder().type(column_type::varchar).build_coldef()});
        st.sequence_number() = 42;

        // Put something in shared_fields, simulating a previous read
        chan.shared_fields().push_back(field_view("prev"));
    }

    test_stream& stream() noexcept { return get_stream(chan); }
};

BOOST_AUTO_TEST_CASE(eof)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(create_eof_frame(42, ok_builder().affected_rows(1).info("1st").build()));

            rows_view rv = fns.read_some_rows_dynamic(fix.chan, fix.st).get();
            BOOST_TEST(rv == makerows(1));
            BOOST_TEST_REQUIRE(fix.st.is_complete());
            BOOST_TEST(fix.st.get_affected_rows() == 1u);
            BOOST_TEST(fix.st.get_info() == "1st");
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
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

            rows_view rv = fns.read_some_rows_dynamic(fix.chan, fix.st).get();
            BOOST_TEST(rv == makerows(1, "abc", "von"));
            BOOST_TEST(fix.st.is_reading_rows());
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
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
                .add_bytes(create_eof_frame(44, ok_builder().affected_rows(1).info("1st").build()));

            rows_view rv = fns.read_some_rows_dynamic(fix.chan, fix.st).get();
            BOOST_TEST(rv == makerows(1, "abc", "von"));
            BOOST_TEST_REQUIRE(fix.st.is_complete());
            BOOST_TEST(fix.st.get_affected_rows() == 1u);
            BOOST_TEST(fix.st.get_info() == "1st");
            BOOST_TEST(fix.chan.shared_sequence_number() == 0u);  // not used
        }
    }
}

// All the other error cases are already tested in read_some_rows_impl. Spotcheck
BOOST_AUTO_TEST_CASE(error)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;

            // invalid row
            fix.stream().add_bytes(create_frame(42, {0x02, 0xff}));

            fns.read_some_rows_dynamic(fix.chan, fix.st)
                .validate_error_exact(client_errc::incomplete_message);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
