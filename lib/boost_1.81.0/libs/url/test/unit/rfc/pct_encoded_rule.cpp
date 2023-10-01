//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/pct_encoded_rule.hpp>

#include <boost/url/grammar/parse.hpp>
#include <boost/url/rfc/pchars.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

class pct_encoded_rule_test
{
public:
    void
    run()
    {
        // javadoc
        {
            result< pct_string_view > rv = grammar::parse( "Program%20Files", pct_encoded_rule( pchars ) );
            (void)rv;
        }
    }
};

TEST_SUITE(
    pct_encoded_rule_test,
    "boost.url.pct_encoded_rule");

} // urls
} // boost
