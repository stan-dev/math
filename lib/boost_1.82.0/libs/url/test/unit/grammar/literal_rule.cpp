//
// Copyright (c) 2022 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/literal_rule.hpp>

#include <boost/url/grammar/parse.hpp>
#include <boost/static_assert.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

BOOST_STATIC_ASSERT(is_rule<literal_rule>::value);

struct literal_rule_test
{
    void
    run()
    {
        // test constexpr
        {
            constexpr auto r = literal_rule("HTTP");
            (void)r;
        }

        // javadoc
        {
            result< string_view > rv = parse( "HTTP", literal_rule( "HTTP" ) );
            (void)rv;
        }

        ok(literal_rule("HTTP"), "HTTP", "HTTP");
        ok(literal_rule("--"), "--", "--");
        bad(literal_rule("http"), "", error::need_more);
        bad(literal_rule("http"), "h", error::need_more);
        bad(literal_rule("http"), "ht", error::need_more);
        bad(literal_rule("http"), "htt", error::need_more);
        bad(literal_rule("http"), "x", error::mismatch);
        bad(literal_rule("http"), "ftp", error::mismatch);
        bad(literal_rule("HTTP"), "http", error::mismatch);
    }
};

TEST_SUITE(
    literal_rule_test,
    "boost.url.grammar.literal_rule");

} // grammar
} // urls
} // boost
