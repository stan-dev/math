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
#include <boost/url/params_ref.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/parse_query.hpp>
#include <boost/url/url.hpp>
#include <boost/static_assert.hpp>
#include <boost/core/ignore_unused.hpp>

#include "test_suite.hpp"

#include <iterator>

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

namespace boost {
namespace urls {

BOOST_STATIC_ASSERT(
    ! std::is_default_constructible<
        params_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_constructible<
        params_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_assignable<
        params_ref>::value);

//------------------------------------------------

#define BIGSTR "123456789012345678901234567890"

struct params_ref_test
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

    // check the sequence
    static
    void
    check(
        params_view const& p,
        std::initializer_list<
            param_pct_view> init)
    {
        if(! BOOST_TEST_EQ(
            p.size(), init.size()))
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

    // check that modification produces
    // the string and correct sequence
    static
    void
    check(
        void(*f)(params_ref),
        string_view s0,
        string_view s1,
        std::initializer_list<
            param_pct_view> init)
    {
        url u;
        {
            auto rv = parse_uri_reference(s0);
            if(! BOOST_TEST(rv.has_value()))
                return;
            u = *rv;
        }
        params_ref ps(u.params());
        f(ps);
        BOOST_TEST_EQ(u.encoded_query(), s1);
        if(! BOOST_TEST_EQ(
                ps.size(), init.size()))
            return;
        check(ps, init);
        {
            auto rv = parse_query(s1);
            if(! BOOST_TEST(rv.has_value()))
                return;
            check(*rv, init);
        }
    }

    static
    void
    check(
        void(*f1)(params_ref),
        void(*f2)(params_ref),
        string_view s0,
        string_view s1,
        std::initializer_list<
            param_pct_view> init)
    {
        check(f1, s0, s1, init);
        check(f2, s0, s1, init);
    }

    //--------------------------------------------

    static
    void
    assign(
        params_ref& p,
        std::initializer_list<
            param_view> init)
    {
        p.assign(
            init.begin(),
            init.end());
    };

    static
    auto
    append(
        params_ref& p,
        std::initializer_list<
            param_view> init) ->
        params_ref::iterator
    {
        return p.append(
            init.begin(),
            init.end());
    };

    static
    auto
    insert(
        params_ref& p,
        params_ref::iterator before,
        std::initializer_list<
            param_view> init) ->
        params_ref::iterator
    {
        return p.insert(
            before,
            init.begin(),
            init.end());
    };

    static
    auto
    replace(
        params_ref& p,
        params_ref::iterator from,
        params_ref::iterator to,
        std::initializer_list<
            param_view> init) ->
        params_ref::iterator
    {
        return p.replace(
            from,
            to,
            init.begin(),
            init.end());
    };

    //--------------------------------------------

    static
    void
    testSpecial()
    {
        // params_ref(params_ref)
        {
            url u("?key=value");
            params_ref p0 = u.params();
            BOOST_TEST_EQ(&p0.url(), &u);
            params_ref p1(p0);
            BOOST_TEST_EQ(&p0.url(), &p1.url());
            check(p0, { {"key","value"} });
            check(p1, { {"key","value"} });
        }

        // params_ref(params_ref, encoding_opts)
        {
            encoding_opts opt;
            opt.space_as_plus = true;
            url u("?key=my+value");
            params_ref p0 = u.params(opt);
            BOOST_TEST_EQ(&p0.url(), &u);
            opt.space_as_plus = false;
            params_ref p1(p0, opt);
            BOOST_TEST_EQ(&p0.url(), &p1.url());
            check(p0, { {"key","my value"} });
            check(p1, { {"key","my+value"} });
        }

        // operator=(params_ref)
        {
            url u0( "?key=value" );
            url u1;
            params_ref p0 = u0.params();
            params_ref p1 = u1.params();
            p1 = p0;
            BOOST_TEST_NE(&p0.url(), &p1.url());
            check(p0, { {"key","value"} });
            check(p1, { {"key","value"} });
        }

        // operator=(initializer_list)
        {
            url u;
            u.params() = {
                {"first", "John"},
                {"last", "Doe"}};
            check(u.params(),{
                {"first", "John"},
                {"last", "Doe"}});
        }

        // operator params_view
        {
            url u;
            params_view qp = u.params();
            BOOST_TEST_EQ(
                qp.buffer().data(), u.buffer().data());
        }
    }

    static
    void
    testObservers()
    {
        // url()
        {
            url u;
            params_ref qp = u.params();
            BOOST_TEST_EQ(&qp.url(), &u);
        }
    }

    static
    void
    testModifiers()
    {
        // clear()
        {
            auto const f = [](params_ref qp)
            {
                qp.clear();
            };
            check(f, "", "", {});
            check(f, "?", "", {});
            check(f, "?first=John&last=Doe", "", {});
        }

        // assign(initializer_list)
        // assign(FwdIt first, FwdIt)
        {
            auto const f = [](params_ref qp)
            {
                qp.assign({ {"first", nullptr}, {"last",""}, {"full", "John Doe"} });
            };
            auto const g = [](params_ref qp)
            {
                assign(qp, { {"first",nullptr}, {"last",""}, {"full", "John Doe"} });
            };
            check(f, g, "", "first&last=&full=John%20Doe",
                { {"first",no_value}, {"last",""}, {"full","John Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                qp.assign({ {BIGSTR, nullptr}, {"last",BIGSTR}, {BIGSTR, BIGSTR} });
            };
            auto const g = [](params_ref qp)
            {
                assign(qp, { {BIGSTR, nullptr}, {"last",BIGSTR}, {BIGSTR, BIGSTR} });
            };
            check(f, g, "", 
                BIGSTR "&last=" BIGSTR "&" BIGSTR "=" BIGSTR,
                { {BIGSTR, nullptr}, {"last",BIGSTR}, {BIGSTR, BIGSTR} });
        }

        // append(param_view)
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.append({"=","&#"});
                BOOST_TEST(is_equal(*it, {"=","&#"}));
            };
            check(f, "?", "&%3D=%26%23", { {"",no_value}, {"=","&#"} });
            check(f, "?key=value", "key=value&%3D=%26%23", { {"key","value"}, {"=","&#"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto it = qp.append({"middle",(*qp.begin()).value});
                BOOST_TEST(is_equal(*it, {"middle","John"}));
            };
            check(f, "?first=John&last=Doe", "first=John&last=Doe&middle=John",
                { {"first","John"}, {"last","Doe"}, {"middle","John"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.append({BIGSTR,BIGSTR});
                BOOST_TEST(is_equal(*it, {BIGSTR,BIGSTR}));
            };
            check(f, "?", "&" BIGSTR "=" BIGSTR, { {"",no_value}, {BIGSTR,BIGSTR} });
            check(f, "?key=value", "key=value&" BIGSTR "=" BIGSTR,
                { {"key","value"}, {BIGSTR,BIGSTR} });
        }

        // append(initializer_list)
        // append(FwdIt first, FwdIt)
        {
            auto const f = [](params_ref qp)
            {
                qp.append({ {"first", nullptr}, {"last",""}, {"full", "John Doe"} });
            };
            auto const g = [](params_ref qp)
            {
                append(qp, { {"first",nullptr}, {"last",""}, {"full", "John Doe"} });
            };
            check(f, g, "", "first&last=&full=John%20Doe",
                { {"first",no_value}, {"last",""}, {"full","John Doe"} });
            check(f, g, "?", "&first&last=&full=John%20Doe",
                { {"",no_value}, {"first",no_value}, {"last",""}, {"full","John Doe"} });
            check(f, g, "?key=value", "key=value&first&last=&full=John%20Doe",
                { {"key","value"}, {"first",no_value}, {"last",""}, {"full","John Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                qp.append({ {BIGSTR, nullptr}, {"last",BIGSTR}, {BIGSTR, BIGSTR} });
            };
            auto const g = [](params_ref qp)
            {
                append(qp, { {BIGSTR, nullptr}, {"last",BIGSTR}, {BIGSTR, BIGSTR} });
            };
            check(f, g, "", BIGSTR "&last=" BIGSTR "&" BIGSTR "=" BIGSTR,
                { {BIGSTR,no_value}, {"last",BIGSTR}, {BIGSTR,BIGSTR} });
        }

        // insert(iterator, param_view)
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto it = qp.insert(std::next(qp.begin(),0),
                    {"middle",(*qp.begin()).value});
                BOOST_TEST(is_equal(*it, {"middle","John"}));
            };
            check(f, "?first=John&last=Doe", "middle=John&first=John&last=Doe",
                { {"middle","John"}, {"first","John"}, {"last","Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto it = qp.insert(std::next(qp.begin(),1),
                    {"middle",(*qp.begin()).value});
                BOOST_TEST(is_equal(*it, {"middle","John"}));
            };
            check(f, "?first=John&last=Doe", "first=John&middle=John&last=Doe",
                { {"first","John"}, {"middle","John"}, {"last","Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto it = qp.insert(std::next(qp.begin(),2),
                    {"middle",(*qp.begin()).value});
                BOOST_TEST(is_equal(*it, {"middle","John"}));
            };
            check(f, "?first=John&last=Doe", "first=John&last=Doe&middle=John",
                { {"first","John"}, {"last","Doe"}, {"middle","John"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),0),
                    {"middle",BIGSTR});
                BOOST_TEST(is_equal(*it, {"middle",BIGSTR}));
            };
            check(f, "?first=John&last=Doe", "middle=" BIGSTR "&first=John&last=Doe",
                { {"middle",BIGSTR}, {"first","John"}, {"last","Doe"} });
        }

        // insert(iterator, initializer_list)
        // insert(iterator, FwdIt first, FwdIt)
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),0),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, qp.begin());
            };
            auto const g = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),0),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, qp.begin());
            };
            check(f, g, "?k1&k2=&k3=v3",
                "first=John&last=Doe&k1&k2=&k3=v3",
                { {"first","John"}, {"last","Doe"}, {"k1",no_value}, {"k2",""}, {"k3","v3"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),1),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, std::next(qp.begin(), 1));
            };
            auto const g = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),1),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, std::next(qp.begin(), 1));
            };
            check(f, g, "?k1&k2=&k3=v3",
                "k1&first=John&last=Doe&k2=&k3=v3",
                { {"k1",no_value}, {"first","John"}, {"last","Doe"}, {"k2",""}, {"k3","v3"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),2),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, std::next(qp.begin(), 2));
            };
            auto const g = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),2),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, std::next(qp.begin(), 2));
            };
            check(f, g, "?k1&k2=&k3=v3",
                "k1&k2=&first=John&last=Doe&k3=v3",
                { {"k1",no_value}, {"k2",""}, {"first","John"}, {"last","Doe"}, {"k3","v3"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),3),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, std::next(qp.begin(), 3));
            };
            auto const g = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),3),
                    { {"first","John"}, {"last","Doe"} });
                BOOST_TEST(is_equal(*it, {"first","John"}));
                BOOST_TEST_EQ(it, std::next(qp.begin(), 3));
            };
            check(f, g, "?k1&k2=&k3=v3",
                "k1&k2=&k3=v3&first=John&last=Doe",
                { {"k1",no_value}, {"k2",""}, {"k3","v3"}, {"first","John"}, {"last","Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),0),
                    { {"first",BIGSTR}, {BIGSTR,"Doe"} });
                BOOST_TEST(is_equal(*it, {"first",BIGSTR}));
                BOOST_TEST_EQ(it, qp.begin());
            };
            auto const g = [](params_ref qp)
            {
                auto it = qp.insert(std::next(qp.begin(),0),
                    { {"first",BIGSTR}, {BIGSTR,"Doe"} });
                BOOST_TEST(is_equal(*it, {"first",BIGSTR}));
                BOOST_TEST_EQ(it, qp.begin());
            };
            check(f, g, "?k1&k2=&k3=v3",
                "first=" BIGSTR "&" BIGSTR "=Doe&k1&k2=&k3=v3",
                { {"first",BIGSTR}, {BIGSTR, "Doe"}, {"k1",no_value}, {"k2",""}, {"k3","v3"} });
        }

        // erase(iterator)
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.erase(std::next(qp.begin(),0));
                BOOST_TEST(is_equal(*it, {"last","Doe"}));
            };
            check(f, "?first=John&last=Doe", "last=Doe", { {"last","Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.erase(std::next(qp.begin(),1));
                BOOST_TEST_EQ(it, qp.end());
            };
            check(f, "?first=John&last=Doe", "first=John", { {"first","John"} });
        }

        // erase(iterator, iterator)
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.erase(
                    std::next(qp.begin(),0),
                    std::next(qp.begin(),2));
                BOOST_TEST(is_equal(*it, {"k2","key"}));
            };
            check(f, "?k0&k1=&k2=key", "k2=key", { {"k2","key"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.erase(
                    std::next(qp.begin(),1),
                    std::next(qp.begin(),3));
                BOOST_TEST_EQ(it, qp.end());
            };
            check(f, "?k0&k1=&k2=key", "k0", { {"k0",no_value} });
        }

        // erase(pct_string_view, ignore_case_param)
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto n = qp.erase((*qp.find_last("k1", ignore_case)).value);
                BOOST_TEST_EQ(n, 2);
            };
            check(f, "?k0&k1=&k2=key&k1=value&k3=4&K1=k1", "k0&k2=key&k3=4&K1=k1",
                { {"k0",no_value}, {"k2","key"}, {"k3","4"}, {"K1","k1"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto n = qp.erase("k1", ignore_case);
                BOOST_TEST_EQ(n, 3);
            };
            check(f, "?k0&k1=&k2=key&k1=value&k3=4&K1=5", "k0&k2=key&k3=4",
                { {"k0",no_value}, {"k2","key"}, {"k3","4"} });
        }

        // replace(iterator, param_view)
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.replace(std::next(qp.begin(),0),
                    {"=","&#"});
                BOOST_TEST(is_equal(*it, {"=","&#"}));
            };
            check(f, "?first=John&last=Doe", "%3D=%26%23&last=Doe",
                { {"=","&#"}, {"last","Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.replace(std::next(qp.begin(),1),
                    {"=","&#"});
                BOOST_TEST(is_equal(*it, {"=","&#"}));
            };
            check(f, "?first=John&last=Doe", "first=John&%3D=%26%23",
                { {"first","John"}, {"=","&#"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto it = qp.replace(std::next(qp.begin(),0),
                    *std::next(qp.begin(),1));
                BOOST_TEST(is_equal(*it, {"last","Doe"}));
            };
            check(f, "?first=John&last=Doe", "last=Doe&last=Doe",
                { {"last","Doe"}, {"last","Doe"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto it = qp.replace(std::next(qp.begin(),1),
                    *std::next(qp.begin(),0));
                BOOST_TEST(is_equal(*it, {"first","John"}));
            };
            check(f, "?first=John&last=Doe", "first=John&first=John",
                { {"first","John"}, {"first","John"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.replace(std::next(qp.begin(),0),
                    {"=",BIGSTR});
                BOOST_TEST(is_equal(*it, {"=",BIGSTR}));
            };
            check(f, "?first=John&last=Doe", "%3D=" BIGSTR "&last=Doe",
                { {"=",BIGSTR}, {"last","Doe"} });
        }

        // replace(iterator, iterator, initializer_list)
        // replace(iterator, iterator, FwdIt, FwdIt)
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.replace(
                    std::next(qp.begin(),0),
                    std::next(qp.begin(),2),
                    {{"=","&#"}});
                BOOST_TEST(is_equal(*it, {"=","&#"}));
            };
            auto const g = [](params_ref qp)
            {
                auto it = replace(qp,
                    std::next(qp.begin(),0),
                    std::next(qp.begin(),2),
                    {{"=","&#"}});
                BOOST_TEST(is_equal(*it, {"=","&#"}));
            };
            check(f, g, "?k0&k1=&k2=key", "%3D=%26%23&k2=key",
                { {"=","&#"}, {"k2","key"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.replace(
                    std::next(qp.begin(),0),
                    std::next(qp.begin(),2),
                    {{"=",BIGSTR}});
                BOOST_TEST(is_equal(*it, {"=",BIGSTR}));
            };
            auto const g = [](params_ref qp)
            {
                auto it = replace(qp,
                    std::next(qp.begin(),0),
                    std::next(qp.begin(),2),
                    {{"=",BIGSTR}});
                BOOST_TEST(is_equal(*it, {"=",BIGSTR}));
            };
            check(f, g, "?k0&k1=&k2=key", "%3D=" BIGSTR "&k2=key",
                { {"=",BIGSTR}, {"k2","key"} });
        }

        // unset(iterator)
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.unset(std::next(qp.begin(),0));
                BOOST_TEST(is_equal(*it, {"k0", no_value}));
            };
            check(f, "?k0&k1=&k2=key", "k0&k1=&k2=key",
                { {"k0",no_value}, {"k1",""}, {"k2","key"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.unset(std::next(qp.begin(),1));
                BOOST_TEST(is_equal(*it, {"k1", no_value}));
            };
            check(f, "?k0&k1=&k2=key", "k0&k1&k2=key",
                { {"k0",no_value}, {"k1",no_value}, {"k2","key"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.unset(std::next(qp.begin(),2));
                BOOST_TEST(is_equal(*it, {"k2", no_value}));
            };
            check(f, "?k0&k1=&k2=key", "k0&k1=&k2",
                { {"k0",no_value}, {"k1",""}, {"k2",no_value} });
        }

        // set(iterator, pct_string_view)
        {
            auto const f = [](params_ref qp)
            {
                // self-intersect
                auto it = qp.set(std::next(qp.begin(),0),
                    (*qp.find("k2")).value);
                BOOST_TEST(is_equal(*it, {"k0", "key"}));
            };
            check(f, "?k0&k1=&k2=key", "k0=key&k1=&k2=key",
                { {"k0","key"}, {"k1",""}, {"k2","key"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.set(std::next(qp.begin(),1), "&#");
                BOOST_TEST(is_equal(*it, {"k1", "&#"}));
            };
            check(f, "?k0&k1=&k2=key", "k0&k1=%26%23&k2=key",
                { {"k0",no_value}, {"k1","&#"}, {"k2","key"} });
        }
        {
            auto const f = [](params_ref qp)
            {
                auto it = qp.set(std::next(qp.begin(),1), BIGSTR);
                BOOST_TEST(is_equal(*it, {"k1", BIGSTR}));
            };
            check(f, "?k0&k1=&k2=key", "k0&k1=" BIGSTR "&k2=key",
                { {"k0",no_value}, {"k1",BIGSTR}, {"k2","key"} });
        }
    }

    static
    void
    testJavadocs()
    {
        // url()
        {
        url u( "?key=value" );

        assert( &u.params().url() == &u );
        }

        // assign(intializer_list)
        {
        url u;

        u.params().assign( {{ "first", "John" }, { "last", "Doe" }} );
        }

        // append(param_view)
        {
        url u;

        u.params().append( { "first", "John" } );
        }

        // append(initializer_list)
        {
        url u;

        u.params().append({ { "first", "John" }, { "last", "Doe" } });
        }

        // erase(iterator)
        {
        url u( "?first=John&last=Doe" );

        params_ref::iterator it = u.params().erase( u.params().begin() );

        assert( u.encoded_query() == "last=Doe" );

        ignore_unused(it);
        }

        // replace(iterator, param_view)
        {
        url u( "?first=John&last=Doe" );

        u.params().replace( u.params().begin(), { "title", "Mr" });

        assert( u.encoded_query() == "title=Mr&last=Doe" );
        }

        // unset(iterator)
        {
        url u( "?first=John&last=Doe" );

        u.params().unset( u.params().begin() );

        assert( u.encoded_query() == "first&last=Doe" );
        }

        // set(iterator, value)
        {
        url u( "?id=42&id=69" );

        u.params().set( u.params().begin(), "none" );

        assert( u.encoded_query() == "id=none&id=69" );
        }

        // set(string_view, string_view)
        {
        url u( "?id=42&id=69" );

        u.params().set( "id", "none" );

        assert( u.params().count( "id" ) == 1 );
        }
    }

    static
    void
    testAll()
    {
        testSpecial();
        testObservers();
        testModifiers();
        testJavadocs();
    }

    void
    run()
    {
        testAll();
    }
};

TEST_SUITE(
    params_ref_test,
    "boost.url.params_ref");

} // urls
} // boost
