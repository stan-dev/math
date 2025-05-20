//
// Copyright (c) 2022 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/delim_rule.hpp>

#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/parse.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct delim_rule_test
{
    void
    run()
    {
        // constexpr
        {
            constexpr auto r = delim_rule('.');
            (void)r;
        }

        // javadoc
        {
            system::result< core::string_view > rv = parse( ".", delim_rule('.') );

            (void)rv;
        }

        ok(delim_rule('$'), "$", "$");
        bad(delim_rule('$'), "", error::need_more);
        bad(delim_rule('$'), "x", error::mismatch);

        ok(delim_rule(alpha_chars), "a", "a");
        bad(delim_rule(alpha_chars), "", error::need_more);
        bad(delim_rule(alpha_chars), "1", error::mismatch);
    }
};

TEST_SUITE(
    delim_rule_test,
    "boost.url.grammar.delim_rule");

} // grammar
} // urls
} // boost
