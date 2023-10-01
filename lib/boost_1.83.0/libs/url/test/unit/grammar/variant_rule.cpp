//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/variant_rule.hpp>

#include <boost/url/grammar/delim_rule.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/rfc/absolute_uri_rule.hpp>
#include <boost/url/rfc/authority_rule.hpp>
#include <boost/url/rfc/origin_form_rule.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct variant_rule_test
{
    void
    run()
    {
        // constexpr
        constexpr auto r =
            variant_rule(
                delim_rule('('),
                delim_rule(')'));

        // javadoc
        {
            system::result< variant2::variant< url_view, url_view, authority_view, core::string_view > > rv = grammar::parse(
                "/index.html?width=full",
                variant_rule(
                    origin_form_rule,
                    absolute_uri_rule,
                    authority_rule,
                    delim_rule('*') ) );

            (void)rv;
        }

        ok(r, "(", variant<core::string_view, core::string_view>(
            variant2::in_place_index_t<0>{}, "("));
        ok(r, ")", variant<core::string_view, core::string_view>(
            variant2::in_place_index_t<1>{}, ")"));
        bad(r, "", error::mismatch);
        bad(r, "[]", error::mismatch);
    }
};

TEST_SUITE(
    variant_rule_test,
    "boost.url.grammar.variant_rule");

} // grammar
} // urls
} // boost
