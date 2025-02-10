//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/parse_path.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

struct parse_path_test
{
    void
    testFunctions()
    {
        auto rv = parse_path("/path/to/file.txt");
        BOOST_TEST(rv.has_value());
    }

    void
    testJavadocs()
    {
    }

    void
    run()
    {
        testFunctions();
        testJavadocs();
    }
};

TEST_SUITE(
    parse_path_test,
    "boost.url.parse_path");

} // urls
} // boost
