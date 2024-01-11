//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/parse.hpp>

#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/grammar/dec_octet_rule.hpp>
#include "test_suite.hpp"

namespace boost {
namespace urls {
namespace grammar {

struct parse_test
{
    void
    testRef()
    {
        BOOST_STATIC_ASSERT(is_rule<
            decltype(ref(dec_octet_rule))>::value);
        BOOST_TEST(parse("255",
            ref(dec_octet_rule)).has_value());
    }

    void
    run()
    {
        testRef();
    }
};

TEST_SUITE(
    parse_test,
    "boost.url.grammar.parse");

} // grammar
} // urls
} // boost
