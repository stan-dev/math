//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/statement.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>
#include <boost/mysql/impl/internal/network_algorithms/close_statement.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_unit/create_channel.hpp"
#include "test_unit/create_statement.hpp"
#include "test_unit/test_stream.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::mysql::detail::channel;

BOOST_AUTO_TEST_SUITE(test_close_statement)

using netfun_maker = netfun_maker_fn<void, channel&, const statement&>;

struct
{
    netfun_maker::signature close_statement;
    const char* name;
} all_fns[] = {
    {netfun_maker::sync_errc(&detail::close_statement_impl),           "sync" },
    {netfun_maker::async_errinfo(&detail::async_close_statement_impl), "async"},
};

constexpr std::uint8_t expected_message[] = {0x05, 0x00, 0x00, 0x00, 0x19, 0x03, 0x00, 0x00, 0x00};

struct fixture
{
    channel chan{create_channel()};
    statement stmt{statement_builder().id(3).num_params(2).build()};

    test_stream& stream() noexcept { return get_stream(chan); }
};

BOOST_AUTO_TEST_CASE(success)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;

            // Call the function
            fns.close_statement(fix.chan, fix.stmt).validate_no_error();

            // Verify the message we sent
            BOOST_MYSQL_ASSERT_BUFFER_EQUALS(fix.stream().bytes_written(), expected_message);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_network)
{
    for (auto fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.stream().set_fail_count(fail_count(0, client_errc::static_row_parsing_error));

            // Call the function
            fns.close_statement(fix.chan, fix.stmt)
                .validate_error_exact(client_errc::static_row_parsing_error);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
