//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/static_url.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/url.hpp>
#include <boost/url/url_view.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/static_assert.hpp>

#include "test_suite.hpp"

#include <iostream>
#include <type_traits>

namespace boost {
namespace urls {

struct static_url_test
{
    BOOST_STATIC_ASSERT(
        std::is_default_constructible<
            static_url<10>>::value);

    BOOST_STATIC_ASSERT(
        std::is_copy_constructible<
            static_url<10>>::value);

    BOOST_STATIC_ASSERT(
        std::is_copy_assignable<
            static_url<10>>::value);

    BOOST_STATIC_ASSERT(
        std::is_convertible<
            static_url<10>, url_view>::value);

    BOOST_STATIC_ASSERT(
        std::is_convertible<
            static_url<10>, url>::value);

    void
    testSpecial()
    {
        // static_url()
        {
            static_url<1024> u;
            BOOST_TEST_EQ(*u.c_str(), '\0');
            BOOST_TEST(u.buffer().empty());
        }

        // static_url(core::string_view)
        {
            // invalid
            BOOST_TEST_THROWS(
                static_url<1024>("$:$"),
                system::system_error);

            // too large
            BOOST_TEST_THROWS(
                static_url<2>("http://www.example.com"),
                system::system_error);

            // URI
            BOOST_TEST_NO_THROW(
                static_url<128>("http://www.example.com"));

            // relative-ref
            BOOST_TEST_NO_THROW(
                static_url<128>("path/to/file.txt"));
        }

        // static_url(static_url)
        // static_url(url_view_base)
        {
            {
                core::string_view s = "/path/to/file.txt";
                static_url<20> u0(s);
                static_url<20> u1(u0);
                BOOST_TEST_EQ(u0.buffer(), u1.buffer());
                BOOST_TEST_NE(u0.buffer().data(), s.data());
                BOOST_TEST_NE(u1.buffer().data(), s.data());
            }
            {
                core::string_view s = "/path/to/file.txt";
                static_url<24> u0(s);
                static_url<20> u1(u0);
                BOOST_TEST_EQ(u0.buffer(), u1.buffer());
                BOOST_TEST_NE(u0.buffer().data(), s.data());
                BOOST_TEST_NE(u1.buffer().data(), s.data());
            }
            {
                core::string_view s = "/path/to/file.txt";
                static_url<24> u0(s);
                static_url<32> u1(u0);
                BOOST_TEST_EQ(u0.buffer(), u1.buffer());
            }
            {
                core::string_view s = "/path/to/file.txt";
                static_url<24> u0(s);
                BOOST_TEST_THROWS(
                    static_url<8>(u0), system::system_error);
            }

            BOOST_TEST_EQ(
                static_url<24>(url_view(
                    "/path/to/file.txt")).buffer(),
                    "/path/to/file.txt");
            BOOST_TEST_EQ(
                static_url<24>(url(
                    "/path/to/file.txt")).buffer(),
                    "/path/to/file.txt");
        }

        // operator=(static_url)
        // operator=(url_view_base)
        {
            {
                core::string_view s = "/path/to/file.txt";
                static_url<24> u0(s);
                static_url<20> u1;
                u1 = u0;
                BOOST_TEST_EQ(u0.buffer(), u1.buffer());
                BOOST_TEST_NE(u0.buffer().data(), s.data());
                BOOST_TEST_NE(u1.buffer().data(), s.data());
            }
            {
                core::string_view s = "/path/to/file.txt";
                static_url<24> u0(s);
                static_url<32> u1;
                u1 = u0;
                BOOST_TEST_EQ(u0.buffer(), u1.buffer());
            }
            {
                core::string_view s = "/path/to/file.txt";
                static_url<24> u0(s);
                static_url<8> u1;
                BOOST_TEST_THROWS(u1 = u0, system::system_error);
            }
            {
                static_url<24> u("http://www.example.com");
                u = url_view("/path/to/file.txt");
                BOOST_TEST_EQ(u.buffer(), "/path/to/file.txt");
            }
            {
                static_url<24> u("http://www.example.com");
                u = url("/path/to/file.txt");
                BOOST_TEST_EQ(u.buffer(), "/path/to/file.txt");
            }
        }
    }

    void
    testUrlBase()
    {
        {
            static_url<64> u("http://example.com");
            BOOST_TEST_EQ(u.buffer(), "http://example.com");
            u.reserve(32);
            BOOST_TEST_EQ(u.capacity(), 65);
        }
        {
            static_url<64> u("http://example.com");
            BOOST_TEST_EQ(u.buffer(), "http://example.com");
            u.clear();
            BOOST_TEST_EQ(u.buffer(), "");
            BOOST_TEST_EQ(u.capacity(), 65);
        }
    }

    void
    testOstream()
    {
        {
            static_url<64> u("http://example.com");
            std::stringstream ss;
            ss << u;
            BOOST_TEST(ss.str() ==
                "http://example.com");
        }
    }

    void
    testJavadocs()
    {
        // {class}
        {
    static_url< 1024 > u( "https://www.example.com" );

    ignore_unused(u);
        }

        // static_url()
        {
        static_url< 1024 > u;

        ignore_unused(u);
        }

        // static_url(core::string_view)
        {
        static_url< 1024 > u( "https://www.example.com" );

        ignore_unused(u);
        }

    }

    void
    run()
    {
        testSpecial();
        testUrlBase();
        testOstream();
        testJavadocs();
    }
};

TEST_SUITE(
    static_url_test,
    "boost.url.static_url");

} // urls
} // boost
