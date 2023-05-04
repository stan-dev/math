//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>

#include <boost/mysql/detail/auth/auth_calculator.hpp>

#include "test_common.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::mysql::client_errc;
using boost::mysql::error_code;
using boost::mysql::string_view;

BOOST_AUTO_TEST_SUITE(test_auth_calculator)

// mysql_native_password
// clang-format off
struct mysql_native_password
{
    auth_calculator calc;
    std::uint8_t challenge_buffer [20] {
        0x79, 0x64, 0x3d, 0x12, 0x1d, 0x71, 0x74, 0x47,
        0x5f, 0x48, 0x3e, 0x3e, 0x0b, 0x62, 0x0a, 0x03,
        0x3d, 0x27, 0x3a, 0x4c
    }; // Values snooped using Wireshark
    std::uint8_t expected_buffer [20] {
        0xf1, 0xb2, 0xfb, 0x1c, 0x8d, 0xe7, 0x5d, 0xb8,
        0xeb, 0xa8, 0x12, 0x6a, 0xd1, 0x0f, 0xe9, 0xb1,
        0x10, 0x50, 0xd4, 0x28
    };
    string_view challenge = makesv(challenge_buffer);
    string_view expected = makesv(expected_buffer);
};
// clang-format on

BOOST_FIXTURE_TEST_CASE(non_empty_password_ssl_false, mysql_native_password)
{
    auto err = calc.calculate("mysql_native_password", "root", challenge, false);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == expected);
    BOOST_TEST(calc.plugin_name() == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_ssl_true, mysql_native_password)
{
    auto err = calc.calculate("mysql_native_password", "root", challenge, true);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == expected);
    BOOST_TEST(calc.plugin_name() == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_ssl_false, mysql_native_password)
{
    auto err = calc.calculate("mysql_native_password", "", challenge, false);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == "");
    BOOST_TEST(calc.plugin_name() == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_ssl_true, mysql_native_password)
{
    auto err = calc.calculate("mysql_native_password", "", challenge, false);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == "");
    BOOST_TEST(calc.plugin_name() == "mysql_native_password");
}

BOOST_FIXTURE_TEST_CASE(bad_challenge_length, mysql_native_password)
{
    BOOST_TEST(
        (calc.calculate("mysql_native_password", "password", "", true)) ==
        make_error_code(client_errc::protocol_value_error)
    );
    BOOST_TEST(
        (calc.calculate("mysql_native_password", "password", "bad_challenge", true)) ==
        make_error_code(client_errc::protocol_value_error)
    );
}

// caching_sha2_password
// clang-format off
struct caching_sha2_password_test
{
    auth_calculator calc;
    std::uint8_t challenge_buffer [20] {
        0x3e, 0x3b, 0x4, 0x55, 0x4, 0x70, 0x16, 0x3a,
        0x4c, 0x15, 0x35, 0x3, 0x15, 0x76, 0x73, 0x22,
        0x46, 0x8, 0x18, 0x1
    }; // Values snooped using the MySQL Python connector
    std::uint8_t expected_buffer [32] {
        0xa1, 0xc1, 0xe1, 0xe9, 0x1b, 0xb6, 0x54, 0x4b,
        0xa7, 0x37, 0x4b, 0x9c, 0x56, 0x6d, 0x69, 0x3e,
        0x6, 0xca, 0x7, 0x2, 0x98, 0xac, 0xd1, 0x6,
        0x18, 0xc6, 0x90, 0x38, 0x9d, 0x88, 0xe1, 0x20
    };
    string_view challenge = makesv(challenge_buffer);
    string_view expected = makesv(expected_buffer);
    string_view cleartext_challenge{"\4"};
};
// clang-format on

BOOST_FIXTURE_TEST_CASE(non_empty_password_challenge_auth_ssl_false, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "hola", challenge, false);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == expected);
    BOOST_TEST(calc.plugin_name() == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_challenge_auth_ssl_true, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "hola", challenge, true);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == expected);
    BOOST_TEST(calc.plugin_name() == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_cleartext_auth_ssl_false, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "hola", cleartext_challenge, false);
    BOOST_TEST(err == make_error_code(client_errc::auth_plugin_requires_ssl));
}

BOOST_FIXTURE_TEST_CASE(non_empty_password_cleartext_auth_ssl_true, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "hola", cleartext_challenge, true);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == std::string("hola") + '\0');
    BOOST_TEST(calc.plugin_name() == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_challenge_auth_ssl_false, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "", challenge, false);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == "");
    BOOST_TEST(calc.plugin_name() == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_challenge_auth_ssl_true, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "", challenge, true);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == "");
    BOOST_TEST(calc.plugin_name() == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_cleartext_auth_ssl_false, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "", cleartext_challenge, false);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == "");
    BOOST_TEST(calc.plugin_name() == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(empty_password_cleartext_auth_ssl_true, caching_sha2_password_test)
{
    auto err = calc.calculate("caching_sha2_password", "", cleartext_challenge, true);
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(calc.response() == "");
    BOOST_TEST(calc.plugin_name() == "caching_sha2_password");
}

BOOST_FIXTURE_TEST_CASE(caching_sha2_bad_challenge_length, caching_sha2_password_test)
{
    BOOST_TEST(
        (calc.calculate("caching_sha2_password", "password", "", true)) ==
        make_error_code(client_errc::protocol_value_error)
    );
    BOOST_TEST(
        (calc.calculate("caching_sha2_password", "password", "bad_challenge", true)) ==
        make_error_code(client_errc::protocol_value_error)
    );
}

// Bad authentication plugin
BOOST_AUTO_TEST_CASE(unknown_auth_plugin)
{
    auth_calculator calc;
    BOOST_TEST(
        (calc.calculate("bad_plugin", "password", "challenge", true)) ==
        make_error_code(client_errc::unknown_auth_plugin)
    );
    BOOST_TEST(
        (calc.calculate("", "password", "challenge", true)) ==
        make_error_code(client_errc::unknown_auth_plugin)
    );
}

BOOST_AUTO_TEST_SUITE_END()  // test_auth_calculator
