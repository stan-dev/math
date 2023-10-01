//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/params_view.hpp>

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
    std::is_default_constructible<
        params_view>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_constructible<
        params_view>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_assignable<
        params_view>::value);

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        params_view::iterator>::value);

struct params_view_test
{
    void
    testMembers()
    {
        // params_view()
        {
            params_view qp;
            BOOST_TEST(qp.empty());
            BOOST_TEST_EQ(qp.buffer(), "");
            BOOST_TEST_EQ(qp.size(), 0);
        }

        // params_view(params_view)
        {
            params_view qp0("first=John&last=Doe");
            params_view qp1(qp0);
            BOOST_TEST_EQ(
                qp0.buffer().data(),
                qp1.buffer().data());
        }

        // params_view(core::string_view)
        {
            try
            {
                core::string_view s = "first=John&last=Doe";
                params_view qp(s);
                BOOST_TEST_PASS();
                BOOST_TEST_EQ(
                    qp.buffer().data(), s.data());
                BOOST_TEST_EQ(qp.buffer(), s);
            }
            catch(std::exception const&)
            {
                BOOST_TEST_FAIL();
            }

            // reserved character
            BOOST_TEST_THROWS(params_view("#"), system::system_error);

            // invalid percent-escape
            BOOST_TEST_THROWS(params_view("%"), system::system_error);
            BOOST_TEST_THROWS(params_view("%F"), system::system_error);
            BOOST_TEST_THROWS(params_view("%FX"), system::system_error);
            BOOST_TEST_THROWS(params_view("%%"), system::system_error);
            BOOST_TEST_THROWS(params_view("FA%"), system::system_error);
        }

        // params_view(core::string_view, encoding_opts)
        {
            try
            {
                encoding_opts opt;
                opt.space_as_plus = true;
                core::string_view s = "name=John+Doe";
                params_view qp(s, opt);
                BOOST_TEST_PASS();
                BOOST_TEST_EQ(
                    qp.buffer().data(), s.data());
                BOOST_TEST_EQ(qp.buffer(), s);
                BOOST_TEST_EQ((*qp.find("name")).value, "John Doe");
            }
            catch(std::exception const&)
            {
                BOOST_TEST_FAIL();
            }
        }

        // operator=(params_view)
        {
            params_view qp0("first=John&last=Doe");
            params_view qp1("key=value");
            qp0 = qp1;
            BOOST_TEST_EQ(
                qp0.buffer().data(),
                qp1.buffer().data());
        }

        // ostream
        {
            params_view qp("first=John&last=Doe");
            std::stringstream ss;
            ss << qp;
            BOOST_TEST_EQ(ss.str(),
                "first=John&last=Doe");
        }
    }

    void
    testJavadocs()
    {
        // {class}
        {
    url_view u( "?first=John&last=Doe" );

    params_view p = u.params();

    ignore_unused(p);
        }
    }

    void
    run()
    {
        testMembers();
        testJavadocs();
    }
};

TEST_SUITE(
    params_view_test,
    "boost.url.params_view");

} // urls
} // boost
