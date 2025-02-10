//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/diagnostics.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/create_diagnostics.hpp"
#include "test_common/printing.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_diagnostics)

BOOST_AUTO_TEST_CASE(operator_equals)
{
    // Both empty
    BOOST_TEST(diagnostics() == diagnostics());

    // Same message
    BOOST_TEST(create_server_diag("abc") == create_server_diag("abc"));

    // Empty vs. message
    BOOST_TEST(!(diagnostics() == create_server_diag("abc")));

    // Different messages
    BOOST_TEST(!(create_server_diag("def") == create_server_diag("abc")));

    // Same message, client vs. server
    BOOST_TEST(!(create_server_diag("abc") == create_client_diag("abc")));
}

BOOST_AUTO_TEST_CASE(operator_not_equals)
{
    BOOST_TEST(!(diagnostics() != diagnostics()));
    BOOST_TEST(!(create_server_diag("abc") != create_server_diag("abc")));
    BOOST_TEST(diagnostics() != create_server_diag("abc"));
    BOOST_TEST(create_server_diag("def") != create_server_diag("abc"));
}

BOOST_AUTO_TEST_CASE(message_accessors)
{
    auto diag = create_server_diag("abc");
    BOOST_TEST(diag.client_message() == "");
    BOOST_TEST(diag.server_message() == "abc");

    diag = create_client_diag("def");
    BOOST_TEST(diag.client_message() == "def");
    BOOST_TEST(diag.server_message() == "");

    diag = diagnostics();
    BOOST_TEST(diag.client_message() == "");
    BOOST_TEST(diag.server_message() == "");
}

BOOST_AUTO_TEST_CASE(clear)
{
    // Clearing server diagnostics
    auto diag_server = create_server_diag("abc");
    diag_server.clear();
    BOOST_TEST(diag_server.client_message() == "");
    BOOST_TEST(diag_server.server_message() == "");

    // Clearing client diagnostics
    auto diag_client = create_client_diag("def");
    diag_client.clear();
    BOOST_TEST(diag_client.client_message() == "");
    BOOST_TEST(diag_client.server_message() == "");

    // Client/server is appropriately reset
    BOOST_TEST(diag_server == diag_client);
}

BOOST_AUTO_TEST_SUITE_END()  // test_diagnostics
