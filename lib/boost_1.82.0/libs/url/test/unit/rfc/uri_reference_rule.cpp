//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/uri_reference_rule.hpp>

#include <boost/url/grammar/parse.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {

class uri_reference_rule_test
{
public:
    void
    run()
    {
        // javadoc
        {
            result< url_view > rv = grammar::parse( "ws://echo.example.com/?name=boost#demo", uri_reference_rule );
            (void)rv;
        }
    }
};

TEST_SUITE(
    uri_reference_rule_test,
    "boost.url.uri_reference_rule");

} // urls
} // boost
