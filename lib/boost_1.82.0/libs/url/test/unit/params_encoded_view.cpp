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
#include <boost/url/params_encoded_view.hpp>

#include <boost/url/parse_query.hpp>
#include <boost/url/url_view.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/static_assert.hpp>
#include <type_traits>

#include "test_suite.hpp"

namespace boost {
namespace urls {

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        params_encoded_view>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_constructible<
        params_encoded_view>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_assignable<
        params_encoded_view>::value);

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        params_encoded_view::iterator>::value);

struct params_encoded_view_test
{
    static
    bool
    is_equal(
        param_pct_view const& p0,
        param_pct_view const& p1)
    {
        return
            p0.key == p1.key &&
            p0.has_value == p1.has_value &&
            (! p0.has_value ||
                p0.value == p1.value);
    }

    static
    void
    check(
        params_encoded_view const& p,
        std::initializer_list<
            param_pct_view> init)
    {
        if(! BOOST_TEST_EQ(p.size(), init.size()))
            return;
        auto it0 = p.begin();
        auto it1 = init.begin();
        auto const end = init.end();
        while(it1 != end)
        {
            BOOST_TEST(is_equal(*it0, *it1));
            auto tmp = it0++;
            BOOST_TEST_EQ(++tmp, it0);
            ++it1;
        }
        // reverse
        if(init.size() > 0)
        {
            it0 = p.end();
            it1 = init.end();
            do
            {
                auto tmp = it0--;
                BOOST_TEST_EQ(--tmp, it0);
                --it1;
                BOOST_TEST(is_equal(*it0, *it1));
            }
            while(it1 != init.begin());
        }
    }

    void
    testMembers()
    {
        // params_encoded_view()
        {
            params_encoded_view qp;
            BOOST_TEST(qp.empty());
            BOOST_TEST_EQ(qp.buffer(), "");
            BOOST_TEST_EQ(qp.size(), 0);
        }

        // params_encoded_view(params_encoded_view)
        {
            params_encoded_view qp0("first=John&last=Doe");
            params_encoded_view qp1(qp0);
            BOOST_TEST_EQ(
                qp0.buffer().data(),
                qp1.buffer().data());
        }

        // params_encoded_view(string_view)
        {
            try
            {
                string_view s = "first=John&last=Doe";
                params_encoded_view qp(s);
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
            BOOST_TEST_THROWS(params_encoded_view("#"), system_error);

            // invalid percent-escape
            BOOST_TEST_THROWS(params_encoded_view("%"), system_error);
            BOOST_TEST_THROWS(params_encoded_view("%F"), system_error);
            BOOST_TEST_THROWS(params_encoded_view("%FX"), system_error);
            BOOST_TEST_THROWS(params_encoded_view("%%"), system_error);
            BOOST_TEST_THROWS(params_encoded_view("FA%"), system_error);
        }

        // operator=(params_encoded_view)
        {
            params_encoded_view qp0("first=John&last=Doe");
            params_encoded_view qp1("key=value");
            qp0 = qp1;
            BOOST_TEST_EQ(
                qp0.buffer().data(),
                qp1.buffer().data());
        }

        // operator params_view()
        {
            params_encoded_view qp0("first=John&last=Doe");
            params_view qp1(qp0);
            BOOST_TEST_EQ(
                qp0.buffer().data(),
                qp1.buffer().data());
        }

        // ostream
        {
            params_encoded_view qp("first=John&last=Doe");
            std::stringstream ss;
            ss << qp;
            BOOST_TEST_EQ(ss.str(),
                "first=John&last=Doe");
        }
    }

    void
    testRange()
    {
        using T = params_encoded_view;

        check( T(""), {} );
        check( T("key"), {{"key", no_value}} );
        check( T("key="), {{"key", ""}} );
        check( T("key=value"), {{"key", "value"}} );
        check( T("key=value&"), {{"key", "value"}, {"", no_value}} );

        check( T(
            "first=John"
            "&last=Doe"
            "&k3="
            "&k4"
            "&"), {
                {"first", "John"},
                {"last", "Doe"},
                {"k3", ""},
                {"k4", no_value},
                {}});
        check( T(""), {} );
        check( T("&"), { {}, {} } );
        check( T("key"), { { "key", no_value } } );
        check( T("key="), { { "key", "" } } );
        check( T("key=value"), { { "key", "value" } } );
        check( T("first=John&last=Doe"), { { "first", "John" }, { "last", "Doe" } } );
        check( T("key=value&"), { { "key", "value" }, {} } );
        check( T("&key=value"), { {}, { "key", "value" } } );
    }

    void
    testJavadocs()
    {
        // {class}
        {
    url_view u( "?first=John&last=Doe" );

    params_encoded_view p = u.encoded_params();

    ignore_unused(p);
        }
    }

    void
    run()
    {
        testMembers();
        testRange();
        testJavadocs();
    }
};

TEST_SUITE(
    params_encoded_view_test,
    "boost.url.params_encoded_view");

} // urls
} // boost
