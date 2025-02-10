//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/absolute_uri_rule.hpp>

#include "test_rule.hpp"

#include <iostream>

namespace boost {
namespace urls {

class absolute_uri_rule_test
{
public:
    void
    run()
    {
        // javadoc
        {
            system::result< url_view > rv = grammar::parse( "http://example.com/index.htm?id=1", absolute_uri_rule );
            (void)rv;
        }

        auto const& t = absolute_uri_rule;

        bad(t, "");
        bad(t, ":");
        bad(t, "http://#");
        bad(t, "http://x.y.z/?a=b&c=d&#");
        bad(t, "http://x.y.z/?a=b&c=d&#frag");
        bad(t, "http://x.y.z/#frag");
        bad(t, "http://%");
        bad(t, "http://?%");

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

        // magnet link
        ok(t, "magnet:?xt=urn:btih:d2474e86c");

        // reg-name might have ipv4 prefix
        ok(t, "http://192.168.0.1.3.a");
    }
};

TEST_SUITE(
    absolute_uri_rule_test,
    "boost.url.absolute_uri_rule");

} // urls
} // boost
