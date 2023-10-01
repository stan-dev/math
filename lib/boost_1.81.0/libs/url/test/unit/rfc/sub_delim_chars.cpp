//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/sub_delim_chars.hpp>

#include <boost/url/grammar/parse.hpp>
#include <boost/url/rfc/pct_encoded_rule.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {

struct sub_delim_chars_test
{
    void
    run()
    {
        // javadoc
        {
            result< pct_string_view > rv = grammar::parse( "Program%20Files", pct_encoded_rule( sub_delim_chars ) );
            (void)rv;
        }

        test_char_set(
            sub_delim_chars,
                "!$&()*+,;=\x27");

    }
};

TEST_SUITE(
    sub_delim_chars_test,
    "boost.url.sub_delim_chars");

} // urls
} // boost
