//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/diagnostics.hpp>

#include <boost/test/unit_test.hpp>

#include "create_diagnostics.hpp"
#include "printing.hpp"

using boost::mysql::diagnostics;
using boost::mysql::test::create_diagnostics;

namespace {

BOOST_AUTO_TEST_SUITE(test_diagnostics)

BOOST_AUTO_TEST_CASE(operator_equals)
{
    BOOST_TEST(diagnostics() == diagnostics());
    BOOST_TEST(create_diagnostics("abc") == create_diagnostics("abc"));
    BOOST_TEST(!(diagnostics() == create_diagnostics("abc")));
    BOOST_TEST(!(create_diagnostics("def") == create_diagnostics("abc")));
}

BOOST_AUTO_TEST_CASE(operator_not_equals)
{
    BOOST_TEST(!(diagnostics() != diagnostics()));
    BOOST_TEST(!(create_diagnostics("abc") != create_diagnostics("abc")));
    BOOST_TEST(diagnostics() != create_diagnostics("abc"));
    BOOST_TEST(create_diagnostics("def") != create_diagnostics("abc"));
}

BOOST_AUTO_TEST_CASE(accessors)
{
    auto diag = create_diagnostics("abc");
    BOOST_TEST(diag.server_message() == "abc");
}

BOOST_AUTO_TEST_SUITE_END()  // test_diagnostics

}  // namespace
