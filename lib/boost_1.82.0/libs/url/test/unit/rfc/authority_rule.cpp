//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/authority_rule.hpp>

#include <boost/url/grammar/parse.hpp>

#include "test_rule.hpp"

#include <type_traits>

namespace boost {
namespace urls {

class authority_rule_test
{
public:
    void
    run()
    {
        // javadoc
        {
            result< authority_view > rv = grammar::parse( "user:pass@example.com:8080", authority_rule );
        }

        auto const& r = authority_rule;

        bad(r, "%");
        bad(r, "host:a");

        ok(r, "");
        ok(r, ":");
        ok(r, "me@you.com");
        ok(r, "user:pass@");
        ok(r, "user:1234");

        {
            auto rv = grammar::parse(
                "x:y@e.com:8080",
                authority_rule);
            BOOST_TEST(rv.has_value());
            auto a = *rv;
            BOOST_TEST(a.host_type() ==
                host_type::name);
            BOOST_TEST(
                a.encoded_host() == "e.com");
            if(BOOST_TEST(a.has_port()))
            {
                BOOST_TEST_EQ(a.port(), "8080");
                BOOST_TEST_EQ(a.port_number(), 8080);
            }
            if(BOOST_TEST(a.has_userinfo()))
            {
                BOOST_TEST_EQ(a.encoded_user(), "x");
                if(BOOST_TEST(a.has_password()))
                    BOOST_TEST_EQ(
                        a.encoded_password(), "y");
            }
        }
    }
};

TEST_SUITE(
    authority_rule_test,
    "boost.url.authority_rule");

} // urls
} // boost
