//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/CPPAlliance/url
//

// Test that header file is self-contained.
#include <boost/url/parse_query.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

struct parse_query_test
{
    void
    testParse()
    {
        result<params_encoded_view> rv;

        rv = parse_query( "key=value" );
    }

    void
    testJavadocs()
    {
    }

    void
    run()
    {
        testParse();
        testJavadocs();
    }
};

TEST_SUITE(
    parse_query_test,
    "boost.url.parse_query");

} // urls
} // boost
