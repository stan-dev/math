//
// Copyright (c) 2022 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/rfc/ipv6_address_rule.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {

struct ipv6_address_rule_test
{
    void
    run()
    {
        // javadoc
        {
            system::result< ipv6_address > rv = grammar::parse( "2001:0db8:85a3:0000:0000:8a2e:0370:7334", ipv6_address_rule );
            (void)rv;
        }
    }
};

TEST_SUITE(
    ipv6_address_rule_test,
    "boost.url.host_ipv6_address_rule");

} // urls
} // boost
