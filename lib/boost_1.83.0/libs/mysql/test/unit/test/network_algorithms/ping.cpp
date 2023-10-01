//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>
#include <boost/mysql/impl/internal/network_algorithms/ping.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_unit/create_channel.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::mysql::detail::channel;

BOOST_AUTO_TEST_SUITE(test_ping)

using netfun_maker = netfun_maker_fn<void, channel&>;

struct
{
    netfun_maker::signature ping;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&detail::ping_impl),           "sync" },
    {netfun_maker::async_errinfo(&detail::async_ping_impl), "async"},
};

struct fixture
{
    channel chan{create_channel()};

    test_stream& stream() noexcept { return get_stream(chan); }
};

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(create_ok_frame(1, ok_builder().build()));

            // Call the function
            fns.ping(fix.chan).validate_no_error();

            // Verify the message we sent
            const std::uint8_t expected_message[] = {0x01, 0x00, 0x00, 0x00, 0x0e};
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), expected_message);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_network)
{
    for (auto fns : all_fns)
    {
        for (int i = 0; i <= 1; ++i)
        {
            BOOST_TEST_CONTEXT(fns.name << " in network transfer " << i)
            {
                fixture fix;
                fix.stream().set_fail_count(fail_count(i, common_server_errc::er_aborting_connection));

                // Call the function
                fns.ping(fix.chan).validate_error_exact(common_server_errc::er_aborting_connection);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(error_response)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().add_bytes(err_builder()
                                       .seqnum(1)
                                       .code(common_server_errc::er_bad_db_error)
                                       .message("my_message")
                                       .build_frame());

            // Call the function
            fns.ping(fix.chan).validate_error_exact(common_server_errc::er_bad_db_error, "my_message");
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
