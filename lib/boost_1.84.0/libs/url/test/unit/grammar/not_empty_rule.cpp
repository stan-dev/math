//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/not_empty_rule.hpp>

#include <boost/url/grammar/digit_chars.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/rfc/pct_encoded_rule.hpp>
#include <boost/url/rfc/unreserved_chars.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct not_empty_rule_test
{
    void
    run()
    {
        // constexpr
        {
            constexpr auto r = not_empty_rule(
                pct_encoded_rule( unreserved_chars ));
            (void)r;
        }

        // javadoc
        {
            system::result< pct_string_view > rv = parse( "Program%20Files",
                not_empty_rule( pct_encoded_rule( unreserved_chars ) ) );

            (void)rv;
        }

        constexpr auto r = not_empty_rule(
            pct_encoded_rule(
                grammar::digit_chars));

        ok(r, "0", "0");
        ok(r, "9", "9");

        bad(r, "", error::mismatch);
        bad(r, "%", error::invalid);
        bad(r, "%f", error::invalid);
        bad(r, "%x", error::invalid);

        {
            core::string_view s("$");
            auto it = s.data();
            auto const end = it + s.size();
            auto rv = parse(it, end, r);
            BOOST_TEST(
                rv.has_error() &&
                rv.error() == error::mismatch);
        }
    }
};

TEST_SUITE(
    not_empty_rule_test,
    "boost.url.grammar.not_empty_rule");

} // grammar
} // urls
} // boost
