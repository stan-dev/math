//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/pct_string_view.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

struct pct_string_view_test
{
    using S = pct_string_view;

    void
    testSpecial()
    {
        // pct_string_view()
        {
            BOOST_TEST(S() == "");
        }

        // pct_string_view(char const*)
        {
            BOOST_TEST(S("") == "");
            BOOST_TEST(S("x") == "x");
            BOOST_TEST(S("%25") == "%25");
        }

        // pct_string_view(char const*, std::size_t)
        {
            BOOST_TEST(S("", 0) == "");
            BOOST_TEST(S("x", 1) == "x");
            BOOST_TEST(S("%25", 3) == "%25");
        }

        // operator core::string_view
        {
            auto const f = [](core::string_view)
            {
            };
            f(S());
        }
    }

    void
    testRelation()
    {
        // ==
        BOOST_TEST(S("x") == S("x"));
        BOOST_TEST(S("x") == "x");
        BOOST_TEST("x" == S("x"));
        BOOST_TEST(S("x") == core::string_view("x"));
        BOOST_TEST(core::string_view("x") == S("x"));
        BOOST_TEST(std::string("x") == S("x"));
        BOOST_TEST(S("x") == std::string("x"));
    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
        BOOST_TEST(S("x") == std::string_view("x"));
        BOOST_TEST(std::string_view("x") == S("x"));
    #endif

        // !=
        BOOST_TEST(S("x") != S("y"));
        BOOST_TEST(S("x") != "y");
        BOOST_TEST("x" != S("y"));
        BOOST_TEST(S("x") != core::string_view("y"));
        BOOST_TEST(core::string_view("x") != S("y"));
        BOOST_TEST(std::string("x") != S("y"));
        BOOST_TEST(S("x") != std::string("y"));
    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
        BOOST_TEST(S("x") != std::string_view("y"));
        BOOST_TEST(std::string_view("x") != S("y"));
    #endif

        // <
        BOOST_TEST(S("x") < S("y"));
        BOOST_TEST(S("x") < "y");
        BOOST_TEST("x" < S("y"));
        BOOST_TEST(S("x") < core::string_view("y"));
        BOOST_TEST(core::string_view("x") < S("y"));
        BOOST_TEST(std::string("x") < S("y"));
        BOOST_TEST(S("x") < std::string("y"));
    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
        BOOST_TEST(S("x") < std::string_view("y"));
        BOOST_TEST(std::string_view("x") < S("y"));
    #endif

        // <=
        BOOST_TEST(S("x") <= S("x"));
        BOOST_TEST(S("x") <= "x");
        BOOST_TEST("x" <= S("x"));
        BOOST_TEST(S("x") <= core::string_view("x"));
        BOOST_TEST(core::string_view("x") <= S("x"));
        BOOST_TEST(std::string("x") <= S("x"));
        BOOST_TEST(S("x") <= std::string("x"));
    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
        BOOST_TEST(S("x") <= std::string_view("x"));
        BOOST_TEST(std::string_view("x") <= S("x"));
    #endif

        // >
        BOOST_TEST(S("y") > S("x"));
        BOOST_TEST(S("y") > "x");
        BOOST_TEST("y" > S("x"));
        BOOST_TEST(S("y") > core::string_view("x"));
        BOOST_TEST(core::string_view("y") > S("x"));
        BOOST_TEST(std::string("y") > S("x"));
        BOOST_TEST(S("y") > std::string("x"));
    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
        BOOST_TEST(S("y") > std::string_view("x"));
        BOOST_TEST(std::string_view("y") > S("x"));
    #endif

        // >=
        BOOST_TEST(S("x") >= S("x"));
        BOOST_TEST(S("x") >= "x");
        BOOST_TEST("x" >= S("x"));
        BOOST_TEST(S("x") >= core::string_view("x"));
        BOOST_TEST(core::string_view("x") >= S("x"));
        BOOST_TEST(std::string("x") >= S("x"));
        BOOST_TEST(S("x") >= std::string("x"));
    #if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
        BOOST_TEST(S("x") >= std::string_view("x"));
        BOOST_TEST(std::string_view("x") >= S("x"));
    #endif

    }

    void
    run()
    {
        testSpecial();
        testRelation();
    }
};

TEST_SUITE(
    pct_string_view_test,
    "boost.url.pct_string_view");

} // urls
} // boost

/*

std::string     query()
pct_string_view encoded_query()
                set_query( core::string_view )
                set_encoded_query( pct_string_view )

1. u.set_query( u.query() )                     // works
2. u.set_query( u.encoded_query() )             // encodes the encoding
3. u.set_encoded_query( u.query() )             // sometimes works, sometimes throws
4. u.set_encoded_query( u.set_encoded_query() ) // works

*/