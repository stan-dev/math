//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/tuple_rule.hpp>

#include <boost/url/grammar/dec_octet_rule.hpp>
#include <boost/url/grammar/delim_rule.hpp>
#include <boost/url/grammar/digit_chars.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>
#include <boost/url/rfc/pct_encoded_rule.hpp>
#include <boost/url/rfc/unreserved_chars.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct tuple_rule_test
{
    template<class R>
    static
    void
    ok( core::string_view s,
        R const& r)
    {
        BOOST_TEST(
            parse(s, r).has_value());
    }

    template<class R>
    static
    void
    bad(core::string_view s,
        R const& r)
    {
        BOOST_TEST(
            parse(s, r).has_error());
    }

    void
    testSequence()
    {
        ok("$", tuple_rule(delim_rule('$')));
        ok("$!", tuple_rule(delim_rule('$'), delim_rule('!')));
        bad("$", tuple_rule(delim_rule('!')));
    }

    void
    testSquelch()
    {
        system::result< std::tuple< pct_string_view, core::string_view > > r1 =
            parse(
                "www.example.com:443",
                tuple_rule(
                    pct_encoded_rule(unreserved_chars + '-' + '.'),
                    squelch( delim_rule( ':' ) ),
                    token_rule( digit_chars ) ) );

        system::result< std::tuple<
                pct_string_view, core::string_view, core::string_view > > r2 =
            parse(
                "www.example.com:443",
                tuple_rule(
                    pct_encoded_rule(unreserved_chars + '-' + '.'),
                    delim_rule( ':' ),
                    token_rule( digit_chars ) ) );

        (void)r1;
        (void)r2;
    }

    void
    run()
    {
        // constexpr
        {
            constexpr auto r1 =
                tuple_rule(
                    delim_rule('.'),
                    delim_rule('.'));
            constexpr auto r2 =
                tuple_rule(
                    squelch( delim_rule('.') ),
                    delim_rule('.'));
            (void)r1;
            (void)r2;
        }

        // javadoc
        {
            system::result< std::tuple< unsigned char, unsigned char, unsigned char, unsigned char > > rv =
                parse( "192.168.0.1", 
                    tuple_rule(
                        dec_octet_rule,
                        squelch( delim_rule('.') ),
                        dec_octet_rule,
                        squelch( delim_rule('.') ),
                        dec_octet_rule,
                        squelch( delim_rule('.') ),
                        dec_octet_rule ) );

            (void)rv;
        }
        testSequence();
        testSquelch();
    }
};

TEST_SUITE(
    tuple_rule_test,
    "boost.url.grammar.tuple_rule");

} // grammar
} // urls
} // boost
