//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/hexdig_chars.hpp>

#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct hexdig_chars_test
{
    void
    run()
    {
        // javadoc
        {
            system::result< core::string_view > rv = parse( "8086FC19", token_rule( hexdig_chars ) );

            (void)rv;
        }

        test_char_set(
            hexdig_chars,
            "0123456789"
            "ABCDEF"
            "abcdef");

        // hexdig_value
        for_each_char(
        [](char c)
        {
            if(hexdig_chars(c))
                BOOST_TEST(
                    hexdig_value(c) >= 0);
            else
                BOOST_TEST(
                    hexdig_value(c) < 0);
        });
    }
};

TEST_SUITE(
    hexdig_chars_test,
    "boost.url.grammar.hexdig_chars");

} // grammar
} // urls
} // boost
