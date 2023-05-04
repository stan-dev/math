//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: httqp://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/params_base.hpp>

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
        params_base>::value);

BOOST_STATIC_ASSERT(
    ! std::is_copy_constructible<
        params_base>::value);

BOOST_STATIC_ASSERT(
    ! std::is_copy_assignable<
        params_base>::value);

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        params_base::iterator>::value);

struct params_base_test
{
    static
    bool
    is_equal(
        param_view const& p0,
        param_view const& p1)
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
        params_view const& p,
        std::initializer_list<
            param_view> init)
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

    static
    void
    check(
        string_view s,
        std::initializer_list<
            param_view> init)
    {
        url_view u(s);
        check(u.params(), init);
    }

    //--------------------------------------------

    void
    testObservers()
    {
        // empty()
        {
            {
                url_view u;
                BOOST_TEST(u.params().empty());
                check(u.params(), {});
            }
            {
                url_view u("?");
                BOOST_TEST(! u.params().empty());
                check(u.params(), {{"", no_value}});
            }
            {
                url_view u("?k=v");
                BOOST_TEST(! u.params().empty());
                check(u.params(), {{"k", "v"}});
            }
        }

        // size()
        {
            {
                url_view u;
                BOOST_TEST_EQ(u.params().size(), 0);
            }
            {
                url_view u("?");
                BOOST_TEST_EQ(u.params().size(), 1);
            }
            {
                url_view u("?k=v&x=y");
                BOOST_TEST_EQ(u.params().size(), 2);
            }
            {
                url_view u("?k0=0&k1=1&k2=&k3&k4=4444#f");
                BOOST_TEST_EQ(
                    u.params().size(), 5);
                check(u.params(), {
                    { "k0", "0" },
                    { "k1", "1" },
                    { "k2", "" },
                    { "k3", no_value },
                    { "k4", "4444" }});
            }
        }

        // begin()
        {
            {
                url_view u;
                BOOST_TEST_EQ(
                    u.params().begin(),
                    u.params().begin());
            }
            {
                url_view u("?");
                BOOST_TEST_NE(
                    u.params().begin(),
                    u.params().end());
            }
        }

        // end()
        {
            {
                url_view u;
                BOOST_TEST_EQ(
                    u.params().end(),
                    u.params().end());
            }
            {
                url_view u("?");
                BOOST_TEST_NE(
                    u.params().end(),
                    u.params().begin());
            }
        }

        {
            url_view u0("?x=1&y=2&x=3&z=4");
            url_view u1("?%78=1&%79=2&%78=3&%7a=4");
            params_view p0 = u0.params();
            params_view p1 = u1.params();

            // contains
            BOOST_TEST(p0.contains("x"));
            BOOST_TEST(p1.contains("x"));
            BOOST_TEST(! p0.contains("X"));
            BOOST_TEST(! p1.contains("X"));
            BOOST_TEST(p0.contains("X", ignore_case));
            BOOST_TEST(p1.contains("X", ignore_case));

            // count
            BOOST_TEST_EQ(p0.count("x"), 2);
            BOOST_TEST_EQ(p0.count("X"), 0);
            BOOST_TEST_EQ(p1.count("x"), 2);
            BOOST_TEST_EQ(p1.count("X"), 0);
            BOOST_TEST_EQ(p0.count("X", ignore_case), 2);
            BOOST_TEST_EQ(p1.count("X", ignore_case), 2);

            // find
            BOOST_TEST_EQ(p0.find("x"), p0.begin());
            BOOST_TEST_EQ(p1.find("x"), p1.begin());
            BOOST_TEST_EQ(p0.find("X", ignore_case), p0.begin());
            BOOST_TEST_EQ(p1.find("X", ignore_case), p1.begin());

            // find(from)
            BOOST_TEST_EQ(
                p0.find(std::next(p0.begin()), "x"),
                std::next(p0.begin(), 2));
            BOOST_TEST_EQ(
                p1.find(std::next(p1.begin()), "x"),
                std::next(p1.begin(), 2));
            BOOST_TEST_EQ(
                p0.find(std::next(p0.begin()),
                    "X", ignore_case),
                std::next(p0.begin(), 2));
            BOOST_TEST_EQ(
                p1.find(std::next(p1.begin()),
                    "X", ignore_case),
                std::next(p1.begin(), 2));
        }

        // (various)
        {
            url_view u(
                "?a=1&%62=2&c=3&c=4"
                "&c=5&d=6&e=7&d=8&f=9#f");
            params_view p = u.params();
            BOOST_TEST_EQ(p.count("a"), 1);
            BOOST_TEST_EQ(p.count("b"), 1);
            BOOST_TEST_EQ(p.count("c"), 3);
            BOOST_TEST_EQ(p.count("d"), 2);
            BOOST_TEST_EQ(p.count("e"), 1);
            BOOST_TEST_EQ(p.count("f"), 1);

            BOOST_TEST_EQ(p.count("g"), 0);
            BOOST_TEST_EQ(p.count("A"), 0);
            BOOST_TEST_EQ(p.count("B"), 0);
            BOOST_TEST_EQ(p.count("C"), 0);
            BOOST_TEST_EQ(p.count("D"), 0);
            BOOST_TEST_EQ(p.count("E"), 0);
            BOOST_TEST_EQ(p.count("F"), 0);
            BOOST_TEST_EQ(p.count("G"), 0);

            BOOST_TEST_EQ(p.count("A", ignore_case), 1);
            BOOST_TEST_EQ(p.count("B", ignore_case), 1);
            BOOST_TEST_EQ(p.count("C", ignore_case), 3);
            BOOST_TEST_EQ(p.count("D", ignore_case), 2);
            BOOST_TEST_EQ(p.count("E", ignore_case), 1);
            BOOST_TEST_EQ(p.count("F", ignore_case), 1);
            BOOST_TEST_EQ(p.count("G", ignore_case), 0);

            BOOST_TEST(p.contains("a"));
            BOOST_TEST(p.contains("b"));
            BOOST_TEST(p.contains("c"));
            BOOST_TEST(p.contains("d"));
            BOOST_TEST(p.contains("e"));
            BOOST_TEST(p.contains("f"));
            BOOST_TEST(! p.contains("g"));

            BOOST_TEST(! p.contains("A"));
            BOOST_TEST(! p.contains("B"));
            BOOST_TEST(! p.contains("C"));
            BOOST_TEST(! p.contains("D"));
            BOOST_TEST(! p.contains("E"));
            BOOST_TEST(! p.contains("F"));
            BOOST_TEST(! p.contains("G"));

            BOOST_TEST(p.contains("A", ignore_case));
            BOOST_TEST(p.contains("B", ignore_case));
            BOOST_TEST(p.contains("C", ignore_case));
            BOOST_TEST(p.contains("D", ignore_case));
            BOOST_TEST(p.contains("E", ignore_case));
            BOOST_TEST(p.contains("F", ignore_case));
            BOOST_TEST(! p.contains("G", ignore_case));
        }
    }

    void
    testIterator()
    {
        using T = params_base::iterator;

        // iterator()
        {
            T t0;
            T t1;
            BOOST_TEST_EQ(t0, t1);
        }

        // operator==()
        {
            url_view u;
            BOOST_TEST_EQ(
                u.encoded_params().begin(),
                u.encoded_params().begin());
        }

        // operator!=()
        {
            url_view u("?");
            BOOST_TEST_NE(
                u.encoded_params().begin(),
                u.encoded_params().end());
        }

        // value_type outlives reference
        {
            params_base::value_type v;
            {
                url_view u("/?a=1&bb=22&ccc=333&dddd=4444#f");
                params_view qp = u.encoded_params();
                params_base::reference r = *qp.begin();
                v = params_base::value_type(r);
            }
            BOOST_TEST_EQ(v.key, "a");
            BOOST_TEST_EQ(v.value, "1");
            BOOST_TEST_EQ(v.has_value, true);
        }
    }

    void
    testRange()
    {
        check( "", {} );
        check( "?", { {} } );
        check( "?&", { {}, {} } );
        check( "?key", { { "key", no_value } } );
        check( "?key=", { { "key", "" } } );
        check( "?key=value", { { "key", "value" } } );
        check( "?first=John&last=Doe", { { "first", "John" }, { "last", "Doe" } } );
        check( "?key=value&", { { "key", "value" }, {} } );
        check( "?&key=value", { {}, { "key", "value" } } );
    }

    void
    testJavadocs()
    {
        // value_type
        {
        params_view::value_type qp( *url_view( "?first=John&last=Doe" ).params().find( "first" ) );
        }

        // buffer()
        {
        assert( url_view( "?first=John&last=Doe" ).params().buffer() == "first=John&last=Doe" );
        }

        // empty()
        {
        assert( ! url_view( "?key=value" ).params().empty() );
        }

        // size()
        {
        assert( url_view( "?key=value").params().size() == 1 );
        }

        // contains()
        {
        assert( url_view( "?first=John&last=Doe" ).params().contains( "first" ) );
        }

        // count()
        {
        assert( url_view( "?first=John&last=Doe" ).params().count( "first" ) == 1 );
        }

        // find()
        {
        assert( (*url_view( "?first=John&last=Doe" ).params().find( "First", ignore_case )).value == "John" );
        }

        // find()
        {
        url_view u( "?First=John&Last=Doe" );

        assert( u.params().find( "first" ) != u.params().find( "first", ignore_case ) );

        ignore_unused(u);
        }

        // find_last()
        {
        assert( (*url_view( "?first=John&last=Doe" ).params().find_last( "last" )).value == "Doe" );
        }

        // find_last()
        {
        url_view u( "?First=John&Last=Doe" );

        assert( u.params().find_last( "last" ) != u.params().find_last( "last", ignore_case ) );

        ignore_unused(u);
        }
    }

    void
    run()
    {
        testObservers();
        testIterator();
        testRange();
        testJavadocs();
    }
};

TEST_SUITE(
    params_base_test,
    "boost.url.params_base");

} // urls
} // boost
