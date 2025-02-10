//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/optional_rule.hpp>

#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/not_empty_rule.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>
#include <boost/optional/optional_io.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct optional_rule_test
{
    void
    run()
    {
        // constexpr
        {
            constexpr auto r =
                optional_rule(token_rule(alpha_chars));
            (void)r;
        }

        // javadoc
        {
            system::result< boost::optional< core::string_view > > rv = parse( "", optional_rule( token_rule( alpha_chars ) ) );

            (void)rv;
        }

        constexpr auto r = optional_rule(
            token_rule(alpha_chars));

        ok(r, "");
        ok(r, "a", core::string_view("a"));

        {
            core::string_view s("$");
            auto it = s.data();
            auto const end = it + s.size();
            auto rv = parse(it, end, r);
            if(BOOST_TEST(!rv.has_error()))
                BOOST_TEST(*rv == boost::none);
        }
    }
};

TEST_SUITE(
    optional_rule_test,
    "boost.url.grammar.optional_rule");

} // grammar
} // urls
} // boost
