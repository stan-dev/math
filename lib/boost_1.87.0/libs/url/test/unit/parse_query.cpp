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
        {
            system::result<params_encoded_view> rv;
            rv = parse_query( "key=value" );
            BOOST_TEST( ! rv.has_error() );
            BOOST_TEST( rv->size() == 1 );
            BOOST_TEST( rv->begin()->key == "key" );
            BOOST_TEST( rv->begin()->value == "value" );
        }

        // issue #757
        {
            auto data = std::string("abc=def&ghi=jkl&mno=pqr");
            core::string_view view( data.data(), data.size() - 2 );
            system::result<params_encoded_view> rv;
            rv = parse_query( view );
            BOOST_TEST( ! rv.has_error() );
            params_encoded_view params = *rv;
            BOOST_TEST( params.size() == 3 );
            auto it = params.begin();
            BOOST_TEST_EQ( it->key, "abc" );
            BOOST_TEST_EQ( it->value, "def" );
            ++it;
            BOOST_TEST_EQ( it->key, "ghi" );
            BOOST_TEST_EQ( it->value, "jkl" );
            ++it;
            BOOST_TEST_EQ( it->key, "mno" );
            BOOST_TEST_EQ( it->value, "p" );
            ++it;
            BOOST_TEST( it == params.end() );
        }
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
