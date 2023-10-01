//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/relative_ref_rule.hpp>

#include <boost/url/grammar/parse.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

class relative_ref_rule_test
{
public:
    void
    run()
    {    
        // javadoc
        {
            result< url_view > rv = grammar::parse( "images/dot.gif?v=hide#a", relative_ref_rule );
            (void)rv;
        }
    }
};

TEST_SUITE(
    relative_ref_rule_test,
    "boost.url.relative_ref_rule");

} // urls
} // boost
