//
// Copyright (c) 2022 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/dec_octet_rule.hpp>

#include <boost/url/grammar/parse.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct dec_octet_rule_test
{
    void
    run()
    {
        // test constexpr
        constexpr auto r = dec_octet_rule;

        // javadoc
        {
            result< unsigned char > rv = parse( "255", dec_octet_rule );

            (void)rv;
        }

        ok(r, "0", 0);
        ok(r, "1", 1);
        ok(r, "9", 9);
        ok(r, "99", 99);
        ok(r, "255", 255);

        bad(r, "", error::mismatch);
        bad(r, "x1", error::mismatch);
        bad(r, "01", error::invalid);
        bad(r, "256", error::invalid);
        bad(r, "260", error::invalid);
        bad(r, "2550", error::invalid);
    }
};

TEST_SUITE(
    dec_octet_rule_test,
    "boost.url.grammar.dec_octet_rule");

} // grammar
} // urls
} // boost
