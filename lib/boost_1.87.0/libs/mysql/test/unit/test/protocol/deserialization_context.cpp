//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/impl/internal/protocol/impl/deserialization_context.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstdint>
#include <limits>

#include "operators.hpp"

using namespace boost::mysql::detail;
using boost::mysql::client_errc;
using boost::mysql::error_code;

BOOST_AUTO_TEST_SUITE(test_deserialization_context)

BOOST_AUTO_TEST_CASE(first_last_size)
{
    const std::array<std::uint8_t, 5> buff{
        {1, 2, 3, 4, 5}
    };
    deserialization_context ctx(buff);

    BOOST_TEST(ctx.first() == buff.data());
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 5u);
}

BOOST_AUTO_TEST_CASE(advance)
{
    const std::array<std::uint8_t, 5> buff{
        {1, 2, 3, 4, 5}
    };
    deserialization_context ctx(buff);

    ctx.advance(1);
    BOOST_TEST(ctx.first() == buff.data() + 1u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 4u);

    ctx.advance(0);
    BOOST_TEST(ctx.first() == buff.data() + 1u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 4u);

    ctx.advance(3);
    BOOST_TEST(ctx.first() == buff.data() + 4u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 1u);

    ctx.advance(1);
    BOOST_TEST(ctx.first() == buff.data() + 5u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 0u);
}

BOOST_AUTO_TEST_CASE(rewind)
{
    const std::array<std::uint8_t, 5> buff{
        {1, 2, 3, 4, 5}
    };
    deserialization_context ctx(buff);

    ctx.advance(4);
    BOOST_TEST(ctx.first() == buff.data() + 4u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 1u);

    ctx.rewind(2);
    BOOST_TEST(ctx.first() == buff.data() + 2u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 3u);

    ctx.advance(1);
    BOOST_TEST(ctx.first() == buff.data() + 3u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 2u);

    ctx.rewind(0);
    BOOST_TEST(ctx.first() == buff.data() + 3u);
    BOOST_TEST(ctx.last() == buff.data() + 5u);
    BOOST_TEST(ctx.size() == 2u);
}

BOOST_AUTO_TEST_CASE(enough_size)
{
    const std::array<std::uint8_t, 5> buff{
        {1, 2, 3, 4, 5}
    };
    deserialization_context ctx(buff);

    ctx.advance(1);
    BOOST_TEST(ctx.enough_size(0));
    BOOST_TEST(ctx.enough_size(1));
    BOOST_TEST(ctx.enough_size(3));
    BOOST_TEST(ctx.enough_size(4));
    BOOST_TEST(!ctx.enough_size(5));
    BOOST_TEST(!ctx.enough_size((std::numeric_limits<std::size_t>::max)()));

    ctx.advance(2);
    BOOST_TEST(ctx.enough_size(2));
    BOOST_TEST(!ctx.enough_size(3));

    ctx.advance(2);
    BOOST_TEST(ctx.enough_size(0));
    BOOST_TEST(!ctx.enough_size(1));
}

BOOST_AUTO_TEST_CASE(get_string)
{
    const std::array<std::uint8_t, 5> buff{
        {0x61, 0x62, 0x63, 0x64, 0x65}
    };
    deserialization_context ctx(buff);

    ctx.advance(1);
    BOOST_TEST(ctx.get_string(0) == "");
    BOOST_TEST(ctx.get_string(1) == "b");
    BOOST_TEST(ctx.get_string(2) == "bc");
    BOOST_TEST(ctx.get_string(4) == "bcde");

    ctx.advance(2);
    BOOST_TEST(ctx.get_string(1) == "d");
    BOOST_TEST(ctx.get_string(2) == "de");

    ctx.advance(2);
    BOOST_TEST(ctx.get_string(0) == "");
}

BOOST_AUTO_TEST_CASE(check_extra_bytes)
{
    const std::array<std::uint8_t, 5> buff{
        {1, 2, 3, 4, 5}
    };
    deserialization_context ctx(buff);

    BOOST_TEST(ctx.check_extra_bytes() == error_code(client_errc::extra_bytes));

    ctx.advance(1);
    BOOST_TEST(ctx.check_extra_bytes() == error_code(client_errc::extra_bytes));

    ctx.advance(3);
    BOOST_TEST(ctx.check_extra_bytes() == error_code(client_errc::extra_bytes));

    ctx.advance(1);
    BOOST_TEST(ctx.check_extra_bytes() == error_code());
}

// Spotcheck: everything works even if an empty span is passed
BOOST_AUTO_TEST_CASE(no_data)
{
    deserialization_context ctx({});
    BOOST_TEST(ctx.first() == nullptr);
    BOOST_TEST(ctx.last() == nullptr);
    BOOST_TEST(ctx.size() == 0u);
    BOOST_TEST(ctx.enough_size(0u));
    BOOST_TEST(!ctx.enough_size(1u));
    BOOST_TEST(ctx.check_extra_bytes() == error_code());
}

// Spotcheck: chain deserialize stops if one of the operations fails
BOOST_AUTO_TEST_CASE(chain_deserialize_error)
{
    // A Deserializable that keeps track of calls and allows setting the return value
    struct mock_deserializable
    {
        deserialize_errc will_return;
        bool called{false};

        mock_deserializable(deserialize_errc e) : will_return(e) {}

        deserialize_errc deserialize(deserialization_context&)
        {
            called = true;
            return will_return;
        }
    };

    // Setup
    const std::array<std::uint8_t, 5> buff{
        {1, 2, 3, 4, 5}
    };
    deserialization_context ctx(buff);
    mock_deserializable values[] = {
        {deserialize_errc::ok},
        {deserialize_errc::incomplete_message},
        {deserialize_errc::ok}
    };

    // Call the function
    auto err = ctx.deserialize(values[0], values[1], values[2]);

    // Validate
    BOOST_TEST(err == deserialize_errc::incomplete_message);
    BOOST_TEST(values[0].called);
    BOOST_TEST(values[1].called);
    BOOST_TEST(!values[2].called);
}

// going from deserialize_errc to error code
BOOST_AUTO_TEST_CASE(to_error_code_)
{
    BOOST_TEST(to_error_code(deserialize_errc::ok) == error_code());
    BOOST_TEST(
        to_error_code(deserialize_errc::incomplete_message) == error_code(client_errc::incomplete_message)
    );
    BOOST_TEST(
        to_error_code(deserialize_errc::protocol_value_error) == error_code(client_errc::protocol_value_error)
    );
    BOOST_TEST(
        to_error_code(deserialize_errc::server_unsupported) == error_code(client_errc::server_unsupported)
    );
}

BOOST_AUTO_TEST_SUITE_END()
