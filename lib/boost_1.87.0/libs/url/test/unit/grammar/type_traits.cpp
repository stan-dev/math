//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/type_traits.hpp>

#include <string>
#include <vector>

#include "test_suite.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct type_traits_test
{
    void
    run()
    {
    }
};

TEST_SUITE(
    type_traits_test,
    "boost.url.grammar.type_traits");

} // grammar
} // urls
} // boost
