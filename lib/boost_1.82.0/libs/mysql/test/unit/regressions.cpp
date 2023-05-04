//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/connection.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/awaitable.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/test/unit_test.hpp>

#include <cstdint>

#include "assert_buffer_equals.hpp"
#include "buffer_concat.hpp"
#include "create_message.hpp"
#include "run_coroutine.hpp"
#include "test_connection.hpp"

#ifdef BOOST_ASIO_HAS_CO_AWAIT

using namespace boost::mysql::test;
using boost::asio::buffer;
using boost::mysql::results;

namespace {

BOOST_AUTO_TEST_SUITE(test_regressions)

// Make sure async_query() and friends don't cause side
// effects in the initiation
BOOST_AUTO_TEST_CASE(side_effects_in_initiation)
{
    boost::asio::io_context ctx;
    test_connection conn;
    results result1, result2;

    // Resultsets will be complete as soon as a message is read
    auto ok_packet_1 = create_ok_packet_message(1, 1, 6, 0, 9, "ab");
    auto ok_packet_2 = create_ok_packet_message(1, 2, 0, 0, 0, "uv");
    conn.stream().add_message(ok_packet_2);
    conn.stream().add_message(ok_packet_1);

    // Launch coroutine and wait for completion
    run_coroutine([&]() -> boost::asio::awaitable<void> {
        // Call both queries but don't wait on them yet, so they don't initiate
        auto aw1 = conn.async_query("Q1", result1, boost::asio::use_awaitable);
        auto aw2 = conn.async_query("Q2", result2, boost::asio::use_awaitable);

        // Run them in reverse order
        co_await std::move(aw2);
        co_await std::move(aw1);
    });

    // Check that we wrote Q2's message first, then Q1's
    auto expected = concat_copy(
        create_message(0, {0x03, 0x51, 0x32}),  // query request Q2
        create_message(0, {0x03, 0x51, 0x31})   // query request Q1
    );
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buffer(conn.stream().bytes_written()), buffer(expected));

    // Check that the results got the right ok_packets
    BOOST_TEST(result1.affected_rows() == 1u);
    BOOST_TEST(result2.affected_rows() == 2u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace

#endif  // BOOST_ASIO_HAS_CO_AWAIT
