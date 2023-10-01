//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/parse.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

struct parse_test
{
    void
    run()
    {
        // issue 497
        {
            auto ru = parse_uri_reference("?~");
            BOOST_TEST_NO_THROW(ru);
            BOOST_TEST(ru->query() == "~");
        }
    }
};

TEST_SUITE(
    parse_test,
    "boost.url.parse");

} // urls
} // boost
