//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/token_rule.hpp>

#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/parse.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct token_rule_test
{
    void
    run()
    {
        // constexpr
        constexpr auto r = token_rule(alpha_chars);

        // javadoc
        {
            result< string_view > rv = parse( "abcdef", token_rule( alpha_chars ) );

            (void)rv;
        }

        ok(r, "a", "a");
        bad(r, "", error::need_more);
        bad(r, "1", error::mismatch);
    }
};

TEST_SUITE(
    token_rule_test,
    "boost.url.grammar.token_rule");

} // grammar
} // urls
} // boost
