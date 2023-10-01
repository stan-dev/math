//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include <boost/url.hpp>

#include "test_suite.hpp"

#include <boost/core/ignore_unused.hpp>
#include <iostream>
#include <list>
#include <string>

namespace boost {
namespace urls {

struct doc_3_urls_test
{
    static
    void
    doc_parsing()
    {
        {
        //[code_containers_1_1
            system::result< url_view > r = parse_uri( "https://www.example.com/path/to/file.txt" );
        //]
            ignore_unused(r);
        }
        {
        //[code_containers_1_2
            url_view u1 = parse_uri_reference( "wss://example.com/quote.cgi?symbol=BOOST&currency=USD" ).value();

            url_view u2( "wss://example.com/quote.cgi?symbol=BOOST&currency=USD" );
        //]
            ignore_unused(u1, u2);
        }
    }

    //[code_container_4_1
    auto segs( core::string_view s ) -> std::list< std::string >
    {
        url_view u( s );
        std::list< std::string > seq;
        for( auto seg : u.encoded_segments() )
            seq.push_back( seg.decode() );
        return seq;
    }
    //]

    //[code_container_5_1
    auto parms( core::string_view s ) -> std::list< param >
    {
        url_view u( s );
        std::list< param > seq;
        for( auto qp : u.params() )
            seq.push_back( qp );
        return seq;
    }
    //]

    void
    path_segments()
    {
        auto check = [](
            core::string_view path,
            std::initializer_list<core::string_view> segs,
            bool path_abs)
        {
            auto r1 = parse_path(path);
            BOOST_TEST(r1.has_value());
            auto ss = r1.value();
            BOOST_TEST_EQ(segs.size(), ss.size());
            BOOST_TEST_EQ(path_abs, r1->is_absolute());
            auto it0 = segs.begin();
            auto it1 = ss.begin();
            while (it0 != segs.end())
            {
                BOOST_TEST_EQ(*it0, *it1);
                ++it0;
                ++it1;
            }
        };

        check("", { }, false);
        check("/", { }, true);
        check("./", { "" }, false);
        check("./usr", { "usr" }, false);
        check("/index.htm", { "index.htm" }, true);
        check("/images/cat-pic.gif", { "images", "cat-pic.gif" }, true);
        check("images/cat-pic.gif", { "images", "cat-pic.gif" }, false);
        check("/fast//query", { "fast", "", "query" }, true);
        check("fast//", { "fast", "", "" }, false);
        check("/./", { "" }, true);
        check(".//", { "", "" }, false);
    }

    void
    run()
    {
        // segs()
        {
            url_view u;
            segs("http://example.com/path/to/file.txt");
        }

        path_segments();
    }
};

TEST_SUITE(
    doc_3_urls_test,
    "boost.url.doc.3_urls");

} // urls
} // boost
