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
#include <boost/url/segments_base.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/parse_path.hpp>
#include <boost/url/url_view.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/static_assert.hpp>

#include "test_suite.hpp"

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

namespace boost {
namespace urls {

BOOST_STATIC_ASSERT(
    ! std::is_default_constructible<
        segments_base>::value);

BOOST_STATIC_ASSERT(
    ! std::is_copy_constructible<
        segments_base>::value);

BOOST_STATIC_ASSERT(
    ! std::is_copy_assignable<
        segments_base>::value);

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        segments_base::iterator>::value);

struct segments_base_test
{
    void
    check(
        core::string_view s,
        std::initializer_list<
            segments_base::reference> match)
    {
        auto rv = parse_uri_reference(s);
        if(! BOOST_TEST(rv.has_value()))
            return;
        segments_base const& ps(
            segments_view(rv->encoded_segments()));
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
                segments_base::reference r0(*it0);
                segments_base::reference r1(*it1);
                BOOST_TEST_EQ(r0, r1);
                BOOST_TEST_EQ(*it0, *it1);
                segments_base::value_type v0(*it0);
                segments_base::value_type v1(*it1);
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
                segments_base::reference r0(*it0);
                segments_base::reference r1(*it1);
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
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1910)
#define Ty(x) static_cast<segments_base const&>(segments_view(x))
#else
//#define Ty(x) static_cast<segments_base const&>(x)
#define Ty(x) static_cast<segments_base const&>(segments_view(x))
#endif
        // is_absolute()
        {
            BOOST_TEST(Ty(
                url_view("/index.htm").segments()).is_absolute());
            BOOST_TEST(! Ty(
                url_view("index.htm").segments()).is_absolute());
            BOOST_TEST(Ty(
                segments_view("/index.htm")).is_absolute());
            BOOST_TEST(! Ty(
                segments_view("index.htm")).is_absolute());
            BOOST_TEST(Ty(
                parse_path("/index.htm").value()).is_absolute());
            BOOST_TEST(! Ty(
                parse_path("index.htm").value()).is_absolute());
        }

        // empty()
        {
            BOOST_TEST(Ty(
                url_view("/").segments()).empty());
            BOOST_TEST(! Ty(
                url_view("x").segments()).empty());
            BOOST_TEST(Ty(
                segments_view("/")).empty());
            BOOST_TEST(! Ty(
                segments_view("x")).empty());
            BOOST_TEST(Ty(
                parse_path("/").value()).empty());
            BOOST_TEST(! Ty(
                parse_path("x").value()).empty());
        }

        // size()
        {
            BOOST_TEST_EQ(Ty(
                url_view("/").segments()).size(), 0);
            BOOST_TEST_EQ(Ty(
                url_view("x").segments()).size(), 1);
            BOOST_TEST_EQ(Ty(
                segments_view("/")).size(), 0);
            BOOST_TEST_EQ(Ty(
                segments_view("x")).size(), 1);
            BOOST_TEST_EQ(Ty(
                parse_path("/").value()).size(), 0);
            BOOST_TEST_EQ(Ty(
                parse_path("x").value()).size(), 1);
        }

        // front()
        {
            BOOST_TEST_EQ(Ty(
                url_view("/x/y").segments()).front(), "x");
            BOOST_TEST_EQ(Ty(
                url_view("x/y").segments()).front(), "x");
            BOOST_TEST_EQ(Ty(
                segments_view("/x/y")).front(), "x");
            BOOST_TEST_EQ(Ty(
                segments_view("x/y")).front(), "x");
            BOOST_TEST_EQ(Ty(
                parse_path("/x/y").value()).front(), "x");
            BOOST_TEST_EQ(Ty(
                parse_path("x/y").value()).front(), "x");
        }

        // back()
        {
            BOOST_TEST_EQ(Ty(
                url_view("/x/y").segments()).back(), "y");
            BOOST_TEST_EQ(Ty(
                url_view("x/y").segments()).back(), "y");
            BOOST_TEST_EQ(Ty(
                segments_view("/x/y")).back(), "y");
            BOOST_TEST_EQ(Ty(
                segments_view("x/y")).back(), "y");
            BOOST_TEST_EQ(Ty(
                parse_path("/x/y").value()).back(), "y");
            BOOST_TEST_EQ(Ty(
                parse_path("x/y").value()).back(), "y");
        }
#undef Ty
    }

    void
    testRange()
    {
    /*  Legend

        '.' %2E     '?' %3F
    */
        check( "", {});
        check( "./", { "" });
        check( ".//", { "", "" });
        check( "/", {});
        check( "/./", { "" });
        check( "/.//", { "", "" });
        check( "/%3F", {"?"});
        check( "%2E/", {".", ""});
        check( "./usr", { "usr" });
        check( "/index.htm", { "index.htm" });
        check( "/images/cat-pic.gif", { "images", "cat-pic.gif" });
        check( "images/cat-pic.gif", { "images", "cat-pic.gif" });
        check( "/fast//query", { "fast", "", "query" });
        check( "fast//",  { "fast", "", "" });
    }

    void
    testJavadoc()
    {
        // value_type
        {
        segments_base::value_type ps( url_view( "/path/to/file.txt" ).segments().back() );

        ignore_unused(ps);
        }

        // buffer()
        {
        assert( url_view( "/path/to/file.txt" ).segments().buffer() == "/path/to/file.txt" );
        }

        // is_absolute()
        {
        assert( url_view( "/path/to/file.txt" ).segments().is_absolute() == true );
        }

        // empty()
        {
        assert( ! url_view( "/index.htm" ).segments().empty() );
        }

        // size()
        {
        assert( url_view( "/path/to/file.txt" ).segments().size() == 3 );
        }

        // front()
        {
        assert( url_view( "/path/to/file.txt" ).segments().front() == "path" );
        }

        // back()
        {
        assert( url_view( "/path/to/file.txt" ).segments().back() == "file.txt" );
        }
    }

    void
    run()
    {
        testObservers();
        testRange();
        testJavadoc();
    }
};

TEST_SUITE(
    segments_base_test,
    "boost.url.segments_base");

} // urls
} // boost
