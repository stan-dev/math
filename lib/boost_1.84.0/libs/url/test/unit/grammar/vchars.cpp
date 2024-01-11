//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/vchars.hpp>

#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/token_rule.hpp>

#include "test_rule.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct vchars_test
{
    void
    run()
    {
        // javadoc
        {
            system::result< core::string_view > rv = parse( "JohnDoe", token_rule( vchars ) );

            (void)rv;
        }
    }
};

TEST_SUITE(
    vchars_test,
    "boost.url.grammar.vchars");

} // grammar
} // urls
} // boost
