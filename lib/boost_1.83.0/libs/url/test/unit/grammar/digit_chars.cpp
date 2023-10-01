//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/digit_chars.hpp>

#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct digit_chars_test
{
    void
    run()
    {
        // javadoc
        {
            system::result< core::string_view > rv = parse( "2022", token_rule( digit_chars ) );

            (void)rv;
        }

        test_char_set(
            digit_chars,
            "0123456789");
    }
};

TEST_SUITE(
    digit_chars_test,
    "boost.url.grammar.digit_chars");

} // grammar
} // urls
} // boost
