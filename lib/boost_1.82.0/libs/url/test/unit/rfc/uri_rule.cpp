//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/uri_rule.hpp>

#include "test_rule.hpp"

#include <iostream>

namespace boost {
namespace urls {

class uri_rule_test
{
public:
    void
    run()
    {
        // javadoc
        {
            result< url_view > rv = grammar::parse( "https://www.example.com/index.htm?id=guest#s1", uri_rule );
            (void)rv;
        }

        auto const& t = uri_rule;

        bad(t, "");
        bad(t, ":");
        bad(t, "http://##");

        ok(t, "http:");
        ok(t, "http:x");
        ok(t, "http:x/");
        ok(t, "http:x/x");
        ok(t, "http:x//");
        ok(t, "http://");
        ok(t, "http://x");
        ok(t, "http://x.y.z");
        ok(t, "http://x.y.z/");
        ok(t, "http://x.y.z/?");
        ok(t, "http://x.y.z/?a");
        ok(t, "http://x.y.z/?a=");
        ok(t, "http://x.y.z/?a=b");
        ok(t, "http://x.y.z/?a=b&c=d");
        ok(t, "http://x.y.z/?a=b&c=d&");
        ok(t, "http://x.y.z/?a=b&c=d&#");
        ok(t, "http://x.y.z/?a=b&c=d&#1");
        ok(t, "http://x.y.z/?a=b&c=d&#12");
        ok(t, "http://x.y.z/?a=b&c=d&#12%23");
        ok(t, "http://x.y.z/?a=b&c=d&#12%23%20");
    }
};

TEST_SUITE(
    uri_rule_test,
    "boost.url.uri_rule");

} // urls
} // boost
