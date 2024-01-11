//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/any_stream_impl.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/buffered_stream.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/ssl/stream.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/netfun_helpers.hpp"
#include "test_common/network_result.hpp"
#include "test_unit/test_stream.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
namespace net = boost::asio;
using boost::mysql::error_code;

BOOST_AUTO_TEST_SUITE(test_any_stream_impl)

struct ssl_fixture
{
    net::ssl::context ctx{net::ssl::context::tls_client};
    any_stream_impl<net::ssl::stream<test_stream>> stream{test_stream(), ctx};

    test_stream& inner_stream() noexcept { return stream.stream().lowest_layer(); }
};

// supports_ssl() and ssl activation/deactivation for SSL streams.
// Due to how this is implemented, testing both at the same time is better
BOOST_FIXTURE_TEST_CASE(ssl_supported, ssl_fixture)
{
    // SSL supported but not active
    BOOST_TEST(stream.supports_ssl());
    BOOST_TEST(!stream.ssl_active());

    // Reset SSL when inactive does nothing
    stream.reset_ssl_active();
    BOOST_TEST(stream.supports_ssl());
    BOOST_TEST(!stream.ssl_active());

    // Set SSL active
    stream.set_ssl_active();
    BOOST_TEST(stream.supports_ssl());
    BOOST_TEST(stream.ssl_active());

    // Reset SSL
    stream.reset_ssl_active();
    BOOST_TEST(stream.supports_ssl());
    BOOST_TEST(!stream.ssl_active());
}

// Same for non-SSL streams
BOOST_AUTO_TEST_CASE(ssl_unsupported)
{
    any_stream_impl<test_stream> stream;

    // SSL not supported
    BOOST_TEST(!stream.supports_ssl());
    BOOST_TEST(!stream.ssl_active());

    // Resetting does nothing
    stream.reset_ssl_active();
    BOOST_TEST(!stream.supports_ssl());
    BOOST_TEST(!stream.ssl_active());
}

// reading/writing when !ssl_active() bypasses the SSL stream
constexpr std::array<std::uint8_t, 3> msg{
    {0x01, 0x02, 0x03}
};

BOOST_FIXTURE_TEST_CASE(write_ssl_disabled, ssl_fixture)
{
    error_code ec;

    std::size_t n = stream.write_some(net::buffer(msg), ec);

    BOOST_TEST(ec == error_code());
    BOOST_TEST(n == 3u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(inner_stream().bytes_written(), msg);
}

// async functions in any_stream use any_completion_handler directly,
// which does not work with netfun makers
BOOST_FIXTURE_TEST_CASE(async_write_ssl_disabled, ssl_fixture)
{
    network_result<std::size_t> result;

    stream.async_write_some(net::buffer(msg), as_network_result<std::size_t>(result, stream.get_executor()));
    run_until_completion(stream.get_executor());

    BOOST_TEST(result.get() == 3u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(inner_stream().bytes_written(), msg);
}

BOOST_FIXTURE_TEST_CASE(read_ssl_disabled, ssl_fixture)
{
    error_code ec;
    std::array<std::uint8_t, 3> read_buff{};
    inner_stream().add_bytes(msg);

    std::size_t n = stream.read_some(net::buffer(read_buff), ec);

    BOOST_TEST(ec == error_code());
    BOOST_TEST(n == 3u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(read_buff, msg);
}

BOOST_FIXTURE_TEST_CASE(async_read_ssl_disabled, ssl_fixture)
{
    network_result<std::size_t> result;
    std::array<std::uint8_t, 3> read_buff{};
    inner_stream().add_bytes(msg);

    stream.async_read_some(
        net::buffer(read_buff),
        as_network_result<std::size_t>(result, stream.get_executor())
    );
    run_until_completion(stream.get_executor());

    BOOST_TEST(result.get() == 3u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(read_buff, msg);
}

BOOST_AUTO_TEST_SUITE_END()
