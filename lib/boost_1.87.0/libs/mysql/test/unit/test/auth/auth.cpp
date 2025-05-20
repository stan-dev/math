//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>

#include <boost/mysql/impl/internal/auth/auth.hpp>

#include "test_common/assert_buffer_equals.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::string_view;

BOOST_AUTO_TEST_SUITE(test_auth)

struct fixture
{
    auth_response resp{
        {0x05, 0x01, 0x02},
        "plugin_not_cleared"
    };
};

BOOST_AUTO_TEST_SUITE(mysql_native_password)

// Values snooped using Wireshark
constexpr std::uint8_t challenge[20] = {
    0x79, 0x64, 0x3d, 0x12, 0x1d, 0x71, 0x74, 0x47, 0x5f, 0x48,
    0x3e, 0x3e, 0x0b, 0x62, 0x0a, 0x03, 0x3d, 0x27, 0x3a, 0x4c,
};
constexpr std::uint8_t expected[20] = {
    0xf1, 0xb2, 0xfb, 0x1c, 0x8d, 0xe7, 0x5d, 0xb8, 0xeb, 0xa8,
    0x12, 0x6a, 0xd1, 0x0f, 0xe9, 0xb1, 0x10, 0x50, 0xd4, 0x28,
};

BOOST_FIXTURE_TEST_CASE(non_empty_password_ssl_false, fixture)
{
    auto err = compute_auth_response("mysql_native_password", "root", challenge, false, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, expected);
    BOOST_TEST(resp.plugin_name == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_ssl_true, fixture)
{
    auto err = compute_auth_response("mysql_native_password", "root", challenge, true, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, expected);
    BOOST_TEST(resp.plugin_name == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_ssl_false, fixture)
{
    auto err = compute_auth_response("mysql_native_password", "", challenge, false, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, std::vector<std::uint8_t>());
    BOOST_TEST(resp.plugin_name == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_ssl_true, fixture)
{
    auto err = compute_auth_response("mysql_native_password", "", challenge, false, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, std::vector<std::uint8_t>());
    BOOST_TEST(resp.plugin_name == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(bad_challenge_length, fixture)
{
    std::uint8_t bad_challenge[] = {0x01, 0x02, 0x03};
    auto err = compute_auth_response("mysql_native_password", "password", bad_challenge, true, resp);
    BOOST_TEST(err == client_errc::protocol_value_error);
}

BOOST_FIXTURE_TEST_CASE(bad_challenge_length_empty, fixture)
{
    auto err = compute_auth_response("mysql_native_password", "password", {}, true, resp);
    BOOST_TEST(err == client_errc::protocol_value_error);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(caching_sha2_password)

// Values snooped using the MySQL Python connector
constexpr std::uint8_t challenge[20] = {
    0x3e, 0x3b, 0x4,  0x55, 0x4,  0x70, 0x16, 0x3a, 0x4c, 0x15,
    0x35, 0x3,  0x15, 0x76, 0x73, 0x22, 0x46, 0x8,  0x18, 0x1,
};

constexpr std::uint8_t expected[32] = {
    0xa1, 0xc1, 0xe1, 0xe9, 0x1b, 0xb6, 0x54, 0x4b, 0xa7, 0x37, 0x4b, 0x9c, 0x56, 0x6d, 0x69, 0x3e,
    0x6,  0xca, 0x7,  0x2,  0x98, 0xac, 0xd1, 0x6,  0x18, 0xc6, 0x90, 0x38, 0x9d, 0x88, 0xe1, 0x20,
};

constexpr std::uint8_t cleartext_challenge[1] = {0x04};

BOOST_FIXTURE_TEST_CASE(non_empty_password_challenge_auth_ssl_false, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "hola", challenge, false, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, expected);
    BOOST_TEST(resp.plugin_name == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_challenge_auth_ssl_true, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "hola", challenge, true, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, expected);
    BOOST_TEST(resp.plugin_name == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_cleartext_auth_ssl_false, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "hola", cleartext_challenge, false, resp);
    BOOST_TEST(err == make_error_code(client_errc::auth_plugin_requires_ssl));
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_cleartext_auth_ssl_true, fixture)
{
    constexpr std::uint8_t expected_local[] = {'h', 'o', 'l', 'a', '\0'};
    auto err = compute_auth_response("caching_sha2_password", "hola", cleartext_challenge, true, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, expected_local);
    BOOST_TEST(resp.plugin_name == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_challenge_auth_ssl_false, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "", challenge, false, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, std::vector<std::uint8_t>());
    BOOST_TEST(resp.plugin_name == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_challenge_auth_ssl_true, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "", challenge, true, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, std::vector<std::uint8_t>());
    BOOST_TEST(resp.plugin_name == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_cleartext_auth_ssl_false, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "", cleartext_challenge, false, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, std::vector<std::uint8_t>());
    BOOST_TEST(resp.plugin_name == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_cleartext_auth_ssl_true, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "", cleartext_challenge, true, resp);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(resp.data, std::vector<std::uint8_t>());
    BOOST_TEST(resp.plugin_name == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(caching_sha2_bad_challenge_length, fixture)
{
    constexpr std::uint8_t bad_challenge[] = {0x00, 0x01, 0x02};
    auto err = compute_auth_response("caching_sha2_password", "password", bad_challenge, true, resp);
    BOOST_TEST(err == client_errc::protocol_value_error);
}

BOOST_FIXTURE_TEST_CASE(caching_sha2_bad_challenge_length_empty, fixture)
{
    auto err = compute_auth_response("caching_sha2_password", "password", {}, true, resp);
    BOOST_TEST(err == client_errc::protocol_value_error);
}
BOOST_AUTO_TEST_SUITE_END()

// Bad authentication plugin
BOOST_AUTO_TEST_CASE(unknown_auth_plugin)
{
    auth_response resp;
    constexpr std::uint8_t challenge[] = {0x00, 0x01, 0x02};
    auto err = compute_auth_response("bad_plugin", "password", challenge, true, resp);
    BOOST_TEST(err == make_error_code(client_errc::unknown_auth_plugin));
}

BOOST_AUTO_TEST_CASE(unknown_auth_plugin_empty)
{
    auth_response resp;
    constexpr std::uint8_t challenge[] = {0x00, 0x01, 0x02};
    auto err = compute_auth_response("", "password", challenge, true, resp);
    BOOST_TEST(err == make_error_code(client_errc::unknown_auth_plugin));
}

BOOST_AUTO_TEST_SUITE_END()
