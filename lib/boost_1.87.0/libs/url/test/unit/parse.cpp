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
            auto r = parse_uri_reference("?~");
            BOOST_TEST_NO_THROW(r.value());
            BOOST_TEST(r->query() == "~");
        }
        // issue 665
        {
            {
                auto r = parse_uri_reference("A:\\");
                BOOST_TEST_THROWS(r.value(), system::system_error);
            }
            {
                auto r = parse_uri_reference("A:\"");
                BOOST_TEST_THROWS(r.value(), system::system_error);
            }
        }
        // reg-name might have ipv4 prefix
        {
            BOOST_TEST_NOT(parse_relative_ref("//0.1.0.1%"));
        }
        // parse docs
        {
            system::result< url_view > r = parse_relative_ref( "//www.boost.org/index.html?field=value#downloads" );
            if ( r.has_value() )
            {
                url_view u = *r;
                assert(u.encoded_path() == "/index.html");
            }
        }
        {
            system::result< url_view > r = parse_uri_reference( "https://www.example.com/path/to/file.txt" );
            if ( r.has_value() )
            {
                url_view u = *r;
                assert(u.encoded_path() == "/path/to/file.txt");
            }
        }
    }
};

TEST_SUITE(
    parse_test,
    "boost.url.parse");

} // urls
} // boost
