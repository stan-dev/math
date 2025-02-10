//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/mysql_collations.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/protocol/serialization.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstdint>
#include <vector>

#include "serialization_test.hpp"
#include "test_common/assert_buffer_equals.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_unit/mock_message.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
namespace collations = boost::mysql::mysql_collations;
using boost::span;
using boost::mysql::date;
using boost::mysql::datetime;
using boost::mysql::error_code;
using boost::mysql::field_view;
using boost::mysql::string_view;

BOOST_AUTO_TEST_SUITE(test_serialization)

// spotcheck: multi-frame messages handled correctly by serialize_top_level
BOOST_AUTO_TEST_CASE(serialize_top_level_multiframe)
{
    constexpr std::size_t frame_size = 8u;
    const std::array<std::uint8_t, 11> payload{
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
    };
    const std::vector<std::uint8_t> expected{80, 81, 82, 83, 85, 8, 0, 0, 42, 1, 2,  3,
                                             4,  5,  6,  7,  8,  3, 0, 0, 43, 9, 10, 11};

    std::vector<std::uint8_t> buff{80, 81, 82, 83, 85};
    auto result = serialize_top_level(mock_message{payload}, buff, 42, 0xffff, frame_size);
    BOOST_TEST(result.err == error_code());
    BOOST_TEST(result.seqnum == 44u);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(buff, expected);
}

// spotcheck: max size correctly propagated
BOOST_AUTO_TEST_CASE(serialize_top_level_error_max_size)
{
    const std::array<std::uint8_t, 11> payload{
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
    };
    std::vector<std::uint8_t> buff;
    auto result = serialize_top_level(mock_message{payload}, buff, 42, 8u);
    BOOST_TEST(result.err == boost::mysql::client_errc::max_buffer_size_exceeded);
    BOOST_TEST(result.seqnum == 0u);
}

BOOST_AUTO_TEST_CASE(quit)
{
    quit_command cmd;
    const std::uint8_t serialized[] = {0x01};
    do_serialize_test(cmd, serialized);
}

BOOST_AUTO_TEST_CASE(ping)
{
    ping_command cmd;
    const std::uint8_t serialized[] = {0x0e};
    do_serialize_test(cmd, serialized);
}

BOOST_AUTO_TEST_CASE(reset_connection)
{
    reset_connection_command cmd;
    const std::uint8_t serialized[] = {0x1f};
    do_serialize_test(cmd, serialized);
}

BOOST_AUTO_TEST_CASE(query)
{
    query_command cmd{"show databases"};
    const std::uint8_t serialized[] =
        {0x03, 0x73, 0x68, 0x6f, 0x77, 0x20, 0x64, 0x61, 0x74, 0x61, 0x62, 0x61, 0x73, 0x65, 0x73};
    do_serialize_test(cmd, serialized);
}

BOOST_AUTO_TEST_CASE(prepare_statement)
{
    prepare_stmt_command cmd{"SELECT * from three_rows_table WHERE id = ?"};
    const std::uint8_t serialized[] = {0x16, 0x53, 0x45, 0x4c, 0x45, 0x43, 0x54, 0x20, 0x2a, 0x20, 0x66,
                                       0x72, 0x6f, 0x6d, 0x20, 0x74, 0x68, 0x72, 0x65, 0x65, 0x5f, 0x72,
                                       0x6f, 0x77, 0x73, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x20, 0x57,
                                       0x48, 0x45, 0x52, 0x45, 0x20, 0x69, 0x64, 0x20, 0x3d, 0x20, 0x3f};
    do_serialize_test(cmd, serialized);
}

BOOST_AUTO_TEST_CASE(execute_statement)
{
    constexpr std::uint8_t blob_buffer[] = {0x70, 0x00, 0x01, 0xff};

    struct
    {
        const char* name;
        std::uint32_t stmt_id;
        std::vector<field_view> params;
        std::vector<std::uint8_t> serialized;
    } test_cases[] = {
        // clang-format off
        {
            "uint64_t",
            1,
            make_fv_vector(std::uint64_t(0xabffffabacadae)),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0x08, 0x80, 0xae, 0xad, 0xac, 0xab, 0xff, 0xff, 0xab, 0x00},
        },
        {
            "int64_t",
            1,
            make_fv_vector(std::int64_t(-0xabffffabacadae)),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0x08, 0x00, 0x52, 0x52, 0x53, 0x54, 0x00, 0x00, 0x54, 0xff}
        },
        {
            "string",
            1,
            make_fv_vector(string_view("test")),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0xfe, 0x00, 0x04, 0x74, 0x65, 0x73, 0x74}
        },
        {
            "blob",
            1,
            make_fv_vector(span<const std::uint8_t>(blob_buffer)),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0xfc, 0x00, 0x04, 0x70, 0x00, 0x01, 0xff}
        },
        {
            "float",
            1,
            make_fv_vector(3.14e20f),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0x04, 0x00, 0x01, 0x2d, 0x88, 0x61}
        },
        {
            "double",
            1,
            make_fv_vector(2.1e214),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0x05, 0x00, 0x56, 0xc0, 0xee, 0xa6, 0x95, 0x30, 0x6f, 0x6c}
        },
        {
            "date",
            1,
            make_fv_vector(date(2010u, 9u, 3u)),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0x0a, 0x00, 0x04, 0xda, 0x07, 0x09, 0x03}
        },
        {
            "datetime",
            1,
            make_fv_vector(datetime(2010u, 9u, 3u, 10u, 30u, 59u, 231800u)),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0x0c, 0x00, 0x0b, 0xda, 0x07, 0x09, 0x03, 0x0a, 0x1e, 0x3b,
            0x78, 0x89, 0x03, 0x00}
        },
        {
            "time",
            1,
            make_fv_vector(maket(230, 30, 59, 231800)),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x01, 0x0b, 0x00, 0x0c, 0x00, 0x09, 0x00, 0x00, 0x00, 0x0e, 0x1e,
            0x3b, 0x78, 0x89, 0x03, 0x00}
        },
        {
            "null",
            1,
            make_fv_vector(nullptr),
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x01, 0x06, 0x00}
        },
        {
            "several_params",
            2,
            make_fv_vector(
                std::uint64_t(0xabffffabacadae),
                std::int64_t(-0xabffffabacadae),
                string_view("test"),
                nullptr,
                2.1e214,
                date(2010u, 9u, 3u),
                datetime(2010u, 9u, 3u, 10u, 30u, 59u, 231800u),
                maket(230, 30, 59, 231800),
                nullptr
            ),
            {0x17, 0x02, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x08, 0x01,
            0x01, 0x08, 0x80, 0x08, 0x00, 0xfe, 0x00, 0x06, 0x00, 0x05, 0x00, 0x0a,
            0x00, 0x0c, 0x00, 0x0b, 0x00, 0x06, 0x00, 0xae, 0xad, 0xac, 0xab, 0xff,
            0xff, 0xab, 0x00, 0x52, 0x52, 0x53, 0x54, 0x00, 0x00, 0x54, 0xff, 0x04,
            0x74, 0x65, 0x73, 0x74, 0x56, 0xc0, 0xee, 0xa6, 0x95, 0x30, 0x6f, 0x6c,
            0x04, 0xda, 0x07, 0x09, 0x03, 0x0b, 0xda, 0x07, 0x09, 0x03, 0x0a, 0x1e,
            0x3b, 0x78, 0x89, 0x03, 0x00, 0x0c, 0x00, 0x09, 0x00, 0x00, 0x00, 0x0e,
            0x1e, 0x3b, 0x78, 0x89, 0x03, 0x00}
        },
        {
            "empty",
            1,
            {},
            {0x17, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00}
        }
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            execute_stmt_command cmd{tc.stmt_id, tc.params};
            do_serialize_test(cmd, tc.serialized);
        }
    }
}

BOOST_AUTO_TEST_CASE(close_statement)
{
    close_stmt_command cmd{1};
    const std::uint8_t serialized[] = {0x19, 0x01, 0x00, 0x00, 0x00};
    do_serialize_test(cmd, serialized);
}

BOOST_AUTO_TEST_CASE(login_request_)
{
    constexpr std::array<std::uint8_t, 20> auth_data{
        {0xfe, 0xc6, 0x2c, 0x9f, 0xab, 0x43, 0x69, 0x46, 0xc5, 0x51,
         0x35, 0xa5, 0xff, 0xdb, 0x3f, 0x48, 0xe6, 0xfc, 0x34, 0xc9}
    };

    constexpr std::uint32_t caps = CLIENT_LONG_PASSWORD | CLIENT_LONG_FLAG | CLIENT_LOCAL_FILES |
                                   CLIENT_PROTOCOL_41 | CLIENT_INTERACTIVE | CLIENT_TRANSACTIONS |
                                   CLIENT_SECURE_CONNECTION | CLIENT_MULTI_STATEMENTS | CLIENT_MULTI_RESULTS |
                                   CLIENT_PS_MULTI_RESULTS | CLIENT_PLUGIN_AUTH | CLIENT_CONNECT_ATTRS |
                                   CLIENT_PLUGIN_AUTH_LENENC_CLIENT_DATA |
                                   CLIENT_CAN_HANDLE_EXPIRED_PASSWORDS | CLIENT_SESSION_TRACK |
                                   CLIENT_DEPRECATE_EOF;

    struct
    {
        const char* name;
        login_request value;
        std::vector<std::uint8_t> serialized;
    } test_cases[] = {
        {
         "without_db", {
                capabilities(caps),
                16777216,  // max packet size
                collations::utf8_general_ci,
                "root",  // username
                auth_data,
                "",                       // database; irrelevant, not using connect with DB capability
                "mysql_native_password",  // auth plugin name
            }, {0x85, 0xa6, 0xff, 0x01, 0x00, 0x00, 0x00, 0x01, 0x21, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
             0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
             0x72, 0x6f, 0x6f, 0x74, 0x00, 0x14, 0xfe, 0xc6, 0x2c, 0x9f, 0xab, 0x43, 0x69, 0x46, 0xc5, 0x51,
             0x35, 0xa5, 0xff, 0xdb, 0x3f, 0x48, 0xe6, 0xfc, 0x34, 0xc9, 0x6d, 0x79, 0x73, 0x71, 0x6c, 0x5f,
             0x6e, 0x61, 0x74, 0x69, 0x76, 0x65, 0x5f, 0x70, 0x61, 0x73, 0x73, 0x77, 0x6f, 0x72, 0x64, 0x00},
         },
        {
         "with_db",            {
                capabilities(caps | CLIENT_CONNECT_WITH_DB),
                16777216,  // max packet size
                collations::utf8_general_ci,
                "root",  // username
                auth_data,
                "database",               // DB name
                "mysql_native_password",  // auth plugin name
            },                                    {0x8d, 0xa6, 0xff, 0x01, 0x00, 0x00, 0x00, 0x01, 0x21, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
             0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
             0x00, 0x00, 0x72, 0x6f, 0x6f, 0x74, 0x00, 0x14, 0xfe, 0xc6, 0x2c, 0x9f, 0xab, 0x43, 0x69,
             0x46, 0xc5, 0x51, 0x35, 0xa5, 0xff, 0xdb, 0x3f, 0x48, 0xe6, 0xfc, 0x34, 0xc9, 0x64, 0x61,
             0x74, 0x61, 0x62, 0x61, 0x73, 0x65, 0x00, 0x6d, 0x79, 0x73, 0x71, 0x6c, 0x5f, 0x6e, 0x61,
             0x74, 0x69, 0x76, 0x65, 0x5f, 0x70, 0x61, 0x73, 0x73, 0x77, 0x6f, 0x72, 0x64, 0x00},
         },
    };

    // TODO: test case with collation > 0xff
    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name) { do_serialize_test(tc.value, tc.serialized); }
    }
}

BOOST_AUTO_TEST_CASE(ssl_request_)
{
    constexpr std::uint32_t caps = CLIENT_LONG_FLAG | CLIENT_LOCAL_FILES | CLIENT_PROTOCOL_41 |
                                   CLIENT_INTERACTIVE | CLIENT_SSL | CLIENT_TRANSACTIONS |
                                   CLIENT_SECURE_CONNECTION | CLIENT_MULTI_STATEMENTS | CLIENT_MULTI_RESULTS |
                                   CLIENT_PS_MULTI_RESULTS | CLIENT_PLUGIN_AUTH | CLIENT_CONNECT_ATTRS |
                                   CLIENT_SESSION_TRACK | (1UL << 29);

    // Data
    ssl_request value{
        capabilities(caps),
        0x1000000,  // max packet size
        collations::utf8mb4_general_ci,
    };

    const std::uint8_t serialized[] = {0x84, 0xae, 0x9f, 0x20, 0x00, 0x00, 0x00, 0x01, 0x2d, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

    do_serialize_test(value, serialized);

    // TODO: test case with collation > 0xff
}

BOOST_AUTO_TEST_CASE(auth_switch_response_)
{
    constexpr std::array<std::uint8_t, 20> auth_data{
        {0xba, 0x55, 0x9c, 0xc5, 0x9c, 0xbf, 0xca, 0x06, 0x91, 0xff,
         0xaa, 0x72, 0x59, 0xfc, 0x53, 0xdf, 0x88, 0x2d, 0xf9, 0xcf}
    };

    auth_switch_response value{auth_data};

    constexpr std::array<std::uint8_t, 20> serialized{
        {0xba, 0x55, 0x9c, 0xc5, 0x9c, 0xbf, 0xca, 0x06, 0x91, 0xff,
         0xaa, 0x72, 0x59, 0xfc, 0x53, 0xdf, 0x88, 0x2d, 0xf9, 0xcf}
    };

    do_serialize_test(value, serialized);
}

BOOST_AUTO_TEST_SUITE_END()
