//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/segments_encoded_base.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/parse_path.hpp>
#include <boost/url/url_view.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/static_assert.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

BOOST_STATIC_ASSERT(
    ! std::is_default_constructible<
        segments_encoded_base>::value);

BOOST_STATIC_ASSERT(
    ! std::is_copy_constructible<
        segments_encoded_base>::value);

BOOST_STATIC_ASSERT(
    ! std::is_copy_assignable<
        segments_encoded_base>::value);

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        segments_encoded_base::iterator>::value);

struct segments_encoded_base_test
{
    void
    check(
        core::string_view s,
        std::initializer_list<
            segments_encoded_base::reference> match)
    {
        auto rv = parse_uri_reference(s);
        if(! BOOST_TEST(rv.has_value()))
            return;
        segments_encoded_base const& ps(
            rv->encoded_segments());
        BOOST_TEST_EQ(ps.buffer().data(), s.data());
        BOOST_TEST_EQ(ps.is_absolute(), s.starts_with('/'));
        BOOST_TEST_EQ(ps.empty(), match.size() == 0);
        if(! BOOST_TEST_EQ(ps.size(), match.size()))
            return;
        if(match.size() > 0 && ! ps.empty())
        {
            BOOST_TEST_EQ(ps.front(), *match.begin());
            BOOST_TEST_EQ(ps.back(), *std::prev(match.end()));
        }
        // forward
        {
            auto it0 = ps.begin();
            auto it1 = match.begin();
            auto const end = ps.end();
            while(it0 != end)
            {
                segments_encoded_base::reference r0(*it0);
                segments_encoded_base::reference r1(*it1);
                BOOST_TEST_EQ(r0, r1);
                BOOST_TEST_EQ(*it0, *it1);
                BOOST_TEST_EQ( // arrow
                    it0->size(), it1->size());
                segments_encoded_base::value_type v0(*it0);
                segments_encoded_base::value_type v1(*it1);
                BOOST_TEST_EQ(v0, *it1);
                BOOST_TEST_EQ(v1, *it1);
                auto prev = it0++;
                BOOST_TEST_NE(prev, it0);
                BOOST_TEST_EQ(++prev, it0);
                ++it1;
                BOOST_TEST_EQ(v0, v1);;
            }
        }
        // reverse
        if(match.size() > 0)
        {
            auto const begin = ps.begin();
            auto it0 = ps.end();
            auto it1 = match.end();
            do
            {
                auto prev = it0--;
                BOOST_TEST_NE(prev, it0);
                BOOST_TEST_EQ(--prev, it0);
                --it1;
                segments_encoded_base::reference r0(*it0);
                segments_encoded_base::reference r1(*it1);
                BOOST_TEST_EQ(*it0, *it1);
                BOOST_TEST_EQ(r0, r1);
            }
            while(it0 != begin);
        }
        // ostream
        {
            std::stringstream ss;
            ss << ps;
            BOOST_TEST_EQ(ss.str(),
                rv->encoded_path());
        }
    }

    //--------------------------------------------

    void
    testObservers()
    {
        // is_absolute()
        {
            BOOST_TEST(static_cast<segments_encoded_base const&>(
                url_view("/index.htm").encoded_segments()).is_absolute());
            BOOST_TEST(! static_cast<segments_encoded_base const&>(
                url_view("index.htm").encoded_segments()).is_absolute());
            BOOST_TEST(static_cast<segments_encoded_base const&>(
                segments_encoded_view("/index.htm")).is_absolute());
            BOOST_TEST(! static_cast<segments_encoded_base const&>(
                segments_encoded_view("index.htm")).is_absolute());
            BOOST_TEST(static_cast<segments_encoded_base const&>(
                parse_path("/index.htm").value()).is_absolute());
            BOOST_TEST(! static_cast<segments_encoded_base const&>(
                parse_path("index.htm").value()).is_absolute());
        }

        // empty()
        {
            BOOST_TEST(static_cast<segments_encoded_base const&>(
                url_view("/").encoded_segments()).empty());
            BOOST_TEST(! static_cast<segments_encoded_base const&>(
                url_view("x").encoded_segments()).empty());
            BOOST_TEST(static_cast<segments_encoded_base const&>(
                segments_encoded_view("/")).empty());
            BOOST_TEST(! static_cast<segments_encoded_base const&>(
                segments_encoded_view("x")).empty());
            BOOST_TEST(static_cast<segments_encoded_base const&>(
                parse_path("/").value()).empty());
            BOOST_TEST(! static_cast<segments_encoded_base const&>(
                parse_path("x").value()).empty());
        }

        // size()
        {
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                url_view("/").encoded_segments()).size(), 0);
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                url_view("x").encoded_segments()).size(), 1);
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                segments_encoded_view("/")).size(), 0);
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                segments_encoded_view("x")).size(), 1);
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                parse_path("/").value()).size(), 0);
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                parse_path("x").value()).size(), 1);
        }

        // front()
        {
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                url_view("/x/y").encoded_segments()).front(), "x");
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                url_view("x/y").encoded_segments()).front(), "x");
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                segments_encoded_view("/x/y")).front(), "x");
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                segments_encoded_view("x/y")).front(), "x");
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                parse_path("/x/y").value()).front(), "x");
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                parse_path("x/y").value()).front(), "x");
        }

        // back()
        {
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                url_view("/x/y").encoded_segments()).back(), "y");
            BOOST_TEST_EQ(static_cast<segments_encoded_base const&>(
                url_view("x/y").encoded_segments()).back(), "y");
            BOOST_TEST_EQ(segments_encoded_view("/x/y").back(), "y");
            BOOST_TEST_EQ(segments_encoded_view("x/y").back(), "y");
            BOOST_TEST_EQ(parse_path("/x/y").value().back(), "y");
            BOOST_TEST_EQ(parse_path("x/y").value().back(), "y");
        }
    }

    void
    testRange()
    {
    /*  Legend

        '#' %23     '?' %3F
        '.' %2E     '[' %5B
        '/' %2F     ']' %5D
    */
        check( "", {});
        check( "./", { "" });
        check( ".//", { "", "" });
        check( "/", {});
        check( "/./", { "" });
        check( "/.//", { "", "" });
        check( "/%3F", {"%3F"});
        check( "%2E/", {"%2E", ""});
        check( "./usr", { "usr" });
        check( "/index.htm", { "index.htm" });
        check( "/images/cat-pic.gif", { "images", "cat-pic.gif" });
        check( "images/cat-pic.gif", { "images", "cat-pic.gif" });
        check( "/fast//query", { "fast", "", "query" });
        check( "fast//",  { "fast", "", "" });
    }

    void
    testJavadocs()
    {
        // value_type
        {
        segments_encoded_base::value_type ps( url_view( "/path/to/file.txt" ).encoded_segments().back() );

        ignore_unused(ps);
        }

        // buffer()
        {
        assert( url_view( "/path/to/file.txt" ).encoded_segments().buffer() == "/path/to/file.txt" );
        }

        // is_absolute()
        {
        assert( url_view( "/path/to/file.txt" ).encoded_segments().is_absolute() == true );
        }

        // empty()
        {
        assert( ! url_view( "/index.htm" ).encoded_segments().empty() );
        }

        // size()
        {
        assert( url_view( "/path/to/file.txt" ).encoded_segments().size() == 3 );
        }

        // front()
        {
        assert( url_view( "/path/to/file.txt" ).encoded_segments().front() == "path" );
        }

        // back()
        {
        assert( url_view( "/path/to/file.txt" ).encoded_segments().back() == "file.txt" );
        }
    }

    void
    run()
    {
        testRange();
        testJavadocs();
    }
};

TEST_SUITE(
    segments_encoded_base_test,
    "boost.url.segments_encoded_base");

} // urls
} // boost
