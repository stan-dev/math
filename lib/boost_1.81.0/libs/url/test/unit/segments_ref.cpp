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
#include <boost/url/segments_ref.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/url.hpp>
#include <boost/static_assert.hpp>
#include <boost/core/ignore_unused.hpp>

#include "test_suite.hpp"

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

namespace boost {
namespace urls {

BOOST_STATIC_ASSERT(
    ! std::is_default_constructible<
        segments_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_constructible<
        segments_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_assignable<
        segments_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        segments_ref::iterator>::value);

//------------------------------------------------

#define BIGSTR "123456789012345678901234567890"

struct segments_ref_test
{
    // check that modification produces
    // the string and correct sequence
    static
    void
    check(
        void(*f)(segments_ref),
        string_view s0,
        string_view s1,
        std::initializer_list<
            string_view> init)
    {
        auto rv = parse_uri_reference(s0);
        if(! BOOST_TEST(rv.has_value()))
            return;
        url u = *rv;
        segments_ref ps(u.segments());
        f(ps);
        BOOST_TEST_EQ(u.encoded_path(), s1);
        if(! BOOST_TEST_EQ(
                ps.size(), init.size()))
            return;
        auto it0 = ps.begin();
        auto it1 = init.begin();
        auto const end = ps.end();
        while(it0 != end)
        {
            BOOST_TEST_EQ(*it0, *it1);
            ++it0;
            ++it1;
        }
    }

    static
    void
    check(
        void(*f1)(segments_ref), void(*f2)(segments_ref),
        string_view s0, string_view s1,
        std::initializer_list<
            string_view> init)
    {
        check(f1, s0, s1, init);
        check(f2, s0, s1, init);
    }

    static
    void
    assign(
        segments_ref& ps,
        std::initializer_list<
            string_view> init)
    {
        ps.assign(init.begin(), init.end());
    };

    static
    auto
    insert(
        segments_ref& ps,
        segments_ref::iterator before,
        std::initializer_list<
            string_view> init) ->
        segments_ref::iterator
    {
        return ps.insert(before,
            init.begin(), init.end());
    }

    static
    auto
    replace(
        segments_ref& ps,
        segments_ref::iterator from,
        segments_ref::iterator to,
        std::initializer_list<
            string_view> init) ->
        segments_ref::iterator
    {
        return ps.replace(from, to,
            init.begin(), init.end());
    }

    //--------------------------------------------

    static
    void
    testSpecial()
    {
        // segments_ref(segments_ref)
        {
            url u("/index.htm");
            segments_ref ps0 = u.segments();
            segments_ref ps(ps0);
            BOOST_TEST_EQ(&ps0.url(), &ps.url());
            BOOST_TEST_EQ(
                ps0.url().buffer().data(),
                ps.url().buffer().data());
        }

        // operator=(segments_ref)
        {
            url u1("/index%2Ehtm");
            url u2("/path/to/file.txt");
            segments_ref ps1 = u1.segments();
            segments_ref ps2 = u2.segments();
            BOOST_TEST_NE(
                ps1.buffer().data(),
                ps2.buffer().data());
            ps1 = ps2;
            BOOST_TEST_EQ(
                u1.encoded_path(),
                u2.encoded_path());
            BOOST_TEST_NE(
                ps1.buffer().data(),
                ps2.buffer().data());
        }

        // operator=(segments_view)
        {
            url u1("/index.htm");
            url_view u2("/path/to/file.txt");
            segments_ref ps1 = u1.segments();
            segments_view ps2 = u2.segments();
            BOOST_TEST_NE(
                ps1.buffer().data(),
                ps2.buffer().data());
            ps1 = ps2;
            BOOST_TEST_EQ(
                u1.encoded_path(),
                u2.encoded_path());
            BOOST_TEST_NE(
                ps1.buffer().data(),
                ps2.buffer().data());
        }

        // operator=(initializer_list)
        {
            url u;
            u.segments() = { "path", "to", "file#" };
            BOOST_TEST_EQ(
                u.encoded_path(), "path/to/file%23");
        }

        // operator segments_view()
        {
            url u;
            u.segments() = { "path", "to", "file#" };
            segments_view ps = u.segments();
            auto it = ps.begin();
            BOOST_TEST_EQ(*it++, "path");
            BOOST_TEST_EQ(*it++, "to");
            BOOST_TEST_EQ(*it++, "file#");
        }
    }

    static
    void
    testObservers()
    {
        // url()
        {
            url u0( "/" );
            url u1( "/" );
            BOOST_TEST_EQ(
                &u0.segments().url(), &u0);
            BOOST_TEST_EQ(
                &u1.segments().url(), &u1);
            BOOST_TEST_NE(
                &u0.segments().url(),
                &u1.segments().url());
        }
    }

    static
    void
    testModifiers()
    {
        //
        // clear()
        //

        {
            auto const f = [](segments_ref ps)
            {
                ps.clear();
            };
            check(f, "", "", {} );
            check(f, "/", "/", {});
            check(f, "/index.htm", "/", {});
            check(f, "index.htm", "", {});
            check(f, "/path/to/file.txt", "/", {});
            check(f, "Program%20Files", "", {});
            check(f, "x://y/", "", {});
        }

        //
        // assign(initializer_list)
        // assign(FwdIt, FwdIt)
        //

        {
            auto const f = [](segments_ref ps)
            {
                ps.assign({ "path", "to", "file.txt?" });
            };
            auto const g = [](segments_ref ps)
            {
                assign(ps, { "path", "to", "file.txt?" });
            };
            check(f, g, "", "path/to/file.txt%3F", {"path", "to", "file.txt?"});
            check(f, g, "/", "/path/to/file.txt%3F", {"path", "to", "file.txt?"});
            check(f, g, "/index.htm", "/path/to/file.txt%3F", {"path", "to", "file.txt?"});
            check(f, g, "index.htm", "path/to/file.txt%3F", {"path", "to", "file.txt?"});
            check(f, g, "/path/to/file.txt", "/path/to/file.txt%3F", {"path", "to", "file.txt?"});
            check(f, g, "Program%20Files", "path/to/file.txt%3F", {"path", "to", "file.txt?"});
        }

        {
            auto const f = [](segments_ref ps)
            {
                ps.assign({ BIGSTR, BIGSTR, BIGSTR "?" });
            };
            auto const g = [](segments_ref ps)
            {
                assign(ps, { BIGSTR, BIGSTR, BIGSTR "?" });
            };
            check(f, g, "Program%20Files",
                BIGSTR "/" BIGSTR "/" BIGSTR "%3F",
                {BIGSTR, BIGSTR, BIGSTR "?"});
        }

        //
        // insert(iterator, string_view)
        //

        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(ps.begin(), "");
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "", "./", {""});
            check(f, "/", "/./", {""});
            check(f, "/index.htm", "/.//index.htm", {"", "index.htm"});
            check(f, "index.htm", ".//index.htm", {"", "index.htm"});
            check(f, "path/to/file.txt", ".//path/to/file.txt", {"", "path", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/.//path/to/file.txt", {"", "path", "to", "file.txt"});
            check(f, "Program%20Files", ".//Program%20Files", {"", "Program Files"});
            check(f, "x:", "./", {""});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(ps.begin(), "my seg");
                BOOST_TEST_EQ(*it, "my seg");
            };
            check(f, "", "my%20seg", {"my seg"});
            check(f, "/", "/my%20seg", {"my seg"});
            check(f, "/index.htm", "/my%20seg/index.htm", {"my seg", "index.htm"});
            check(f, "index.htm", "my%20seg/index.htm", {"my seg", "index.htm"});
            check(f, "path/to/file.txt", "my%20seg/path/to/file.txt", {"my seg", "path", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/my%20seg/path/to/file.txt", {"my seg", "path", "to", "file.txt"});
            check(f, "Program%20Files", "my%20seg/Program%20Files", {"my seg", "Program Files"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(std::next(ps.begin(), 1), "my seg");
                BOOST_TEST_EQ(*it, "my seg");
            };
            check(f, "path/to/file.txt", "path/my%20seg/to/file.txt", {"path", "my seg", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/path/my%20seg/to/file.txt", {"path", "my seg", "to", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(ps.end(), "my seg");
                BOOST_TEST_EQ(*it, "my seg");
            };
            check(f, "", "my%20seg", {"my seg"});
            check(f, "/", "/my%20seg", {"my seg"});
            check(f, "/index.htm", "/index.htm/my%20seg", {"index.htm", "my seg"});
            check(f, "index.htm", "index.htm/my%20seg", {"index.htm", "my seg"});
            check(f, "path/to/file.txt", "path/to/file.txt/my%20seg", {"path", "to", "file.txt", "my seg"});
            check(f, "/path/to/file.txt", "/path/to/file.txt/my%20seg", {"path", "to", "file.txt", "my seg"});
            check(f, "Program%20Files", "Program%20Files/my%20seg", {"Program Files", "my seg"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(ps.end(), "");
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "", "./", {""});
            check(f, "/", "/./", {""});
            check(f, "/index.htm", "/index.htm/", {"index.htm", ""});
            check(f, "index.htm", "index.htm/", {"index.htm", ""});
            check(f, "path/to/file.txt", "path/to/file.txt/", {"path", "to", "file.txt", ""});
            check(f, "/path/to/file.txt", "/path/to/file.txt/", {"path", "to", "file.txt", ""});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(std::next(ps.begin(), 1), BIGSTR);
                BOOST_TEST_EQ(*it, BIGSTR);
            };
            check(f, "/path/to/file.txt",
                "/path/" BIGSTR "/to/file.txt",
                {"path", BIGSTR, "to", "file.txt"});
        }

        //
        // insert(iterator, initializer_list)
        // insert(iterator, FwdIt, FwdIt)
        //

        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(ps.begin(), { "u", "v#" });
                BOOST_TEST_EQ(*it, "u");
            };
            auto const g = [](segments_ref ps)
            {
                auto it = insert(ps, ps.begin(), { "u", "v#" });
                BOOST_TEST_EQ(*it, "u");
            };
            check(f, g, "", "u/v%23", {"u", "v#"});
            check(f, g, "/", "/u/v%23", {"u", "v#"});
            check(f, g, "/index.htm", "/u/v%23/index.htm", {"u", "v#", "index.htm"});
            check(f, g, "index.htm", "u/v%23/index.htm", {"u", "v#", "index.htm"});
            check(f, g, "path/to/file.txt", "u/v%23/path/to/file.txt", {"u", "v#", "path", "to", "file.txt"});
            check(f, g, "/path/to/file.txt", "/u/v%23/path/to/file.txt", {"u", "v#", "path", "to", "file.txt"});
            check(f, g, "Program%20Files", "u/v%23/Program%20Files", {"u", "v#", "Program Files"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(ps.begin(), { "", "" });
                BOOST_TEST_EQ(*it, "");
            };
            auto const g = [](segments_ref ps)
            {
                auto it = insert(ps, ps.begin(), { "", "" });
                BOOST_TEST_EQ(*it, "");
            };
            check(f, g, "", ".//", {"", ""});
            check(f, g, "/", "/.//", {"", ""});
            check(f, g, "/index.htm", "/.///index.htm", {"", "", "index.htm"});
            check(f, g, "index.htm", ".///index.htm", {"", "", "index.htm"});
            check(f, g, "path/to/file.txt", ".///path/to/file.txt", {"", "", "path", "to", "file.txt"});
            check(f, g, "/path/to/file.txt", "/.///path/to/file.txt", {"", "", "path", "to", "file.txt"});
            check(f, g, "x", ".///x", {"", "", "x"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.insert(ps.begin(), { BIGSTR, BIGSTR });
                BOOST_TEST_EQ(*it, BIGSTR);
            };
            auto const g = [](segments_ref ps)
            {
                auto it = insert(ps, ps.begin(), { BIGSTR, BIGSTR });
                BOOST_TEST_EQ(*it, BIGSTR);
            };
            check(f, g, "index.htm",
                BIGSTR "/" BIGSTR "/index.htm",
                {BIGSTR, BIGSTR, "index.htm"});
        }

        //
        // erase(iterator)
        //

        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.erase(std::next(ps.begin(), 0));
                if(it != ps.end())
                    BOOST_TEST_EQ(*it, ps.front());
            };
            check(f, "path/to/file.txt", "to/file.txt", {"to", "file.txt"});
            check(f, "/path/to/file.txt", "/to/file.txt", {"to", "file.txt"});
            check(f, "//x/y/", "/./", {""});
            check(f, "/x/", "/./", {""});
            check(f, "x/", "./", {""});
            check(f, "x:.//", "./", {""});
            check(f, ".//:", "./:", {":"});
            check(f, "x:.//:", ":", {":"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.erase(std::next(ps.begin(), 1));
                BOOST_TEST_EQ(*it, "file.txt");
            };
            check(f, "path/to/file.txt", "path/file.txt", {"path", "file.txt"});
            check(f, "/path/to/file.txt", "/path/file.txt", {"path", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.erase(std::next(ps.begin(), 2));
                BOOST_TEST_EQ(it, ps.end());
            };
            check(f, "path/to/file.txt", "path/to", {"path", "to"});
            check(f, "/path/to/file.txt", "/path/to", {"path", "to"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                ps.erase(std::next(ps.begin(), 1));
            };
            check(f, "x://y///", "//", {"", ""});
            check(f, ".///", ".//", {"", ""});
        }

        //
        // erase(iterator, iterator)
        //

        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.erase(
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2));
                BOOST_TEST_EQ(*it, "the");
            };
            check(f, "path/to/the/file.txt", "the/file.txt", {"the", "file.txt"});
            check(f, "/path/to/the/file.txt", "/the/file.txt", {"the", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.erase(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3));
                BOOST_TEST_EQ(*it, "file.txt");
            };
            check(f, "path/to/the/file.txt", "path/file.txt", {"path", "file.txt"});
            check(f, "/path/to/the/file.txt", "/path/file.txt", {"path", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.erase(
                    std::next(ps.begin(), 2),
                    std::next(ps.begin(), 4));
                BOOST_TEST_EQ(it, ps.end());
            };
            check(f, "path/to/the/file.txt", "path/to", {"path", "to"});
            check(f, "/path/to/the/file.txt", "/path/to", {"path", "to"});
        }

        //
        // replace(iterator, string_view)
        //

        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 0), "");
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "path/to/file.txt", ".//to/file.txt", {"", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/.//to/file.txt", {"", "to", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 1), "");
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "path/to/file.txt", "path//file.txt", {"path", "", "file.txt"});
            check(f, "/path/to/file.txt", "/path//file.txt", {"path", "", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 0), "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/file.txt", "test/to/file.txt", {"test", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/test/to/file.txt", {"test", "to", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 1), "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/file.txt", "path/test/file.txt", {"path", "test", "file.txt"});
            check(f, "/path/to/file.txt", "/path/test/file.txt", {"path", "test", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 2), "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/file.txt", "path/to/test", {"path", "to", "test"});
            check(f, "/path/to/file.txt", "/path/to/test", {"path", "to", "test"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(ps.begin(), ":");
                BOOST_TEST_EQ(*it, ":");
            };
            check(f, "path/to/file.txt", "%3A/to/file.txt", {":", "to", "file.txt"});
            check(f, "/:/to/file.txt", "/:/to/file.txt", {":", "to", "file.txt"});
            check(f, "x:path/to/file.txt", ":/to/file.txt", {":", "to", "file.txt"});
            check(f, "x:/path/to/file.txt", "/:/to/file.txt", {":", "to", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 2), BIGSTR);
                BOOST_TEST_EQ(*it, BIGSTR);
            };
            check(f, "path/to/file.txt", "path/to/" BIGSTR, {"path", "to", BIGSTR});
            check(f, "/path/to/file.txt", "/path/to/" BIGSTR, {"path", "to", BIGSTR});
        }

        //
        // replace(iterator, string_view)
        //

        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    string_view(""));
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "path/to/the/file.txt", ".//the/file.txt", {"", "the", "file.txt"});
            check(f, "/path/to/the/file.txt", "/.//the/file.txt", {"", "the", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    "");
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "path/to/the/file.txt", "path//file.txt", {"path", "", "file.txt"});
            check(f, "/path/to/the/file.txt", "/path//file.txt", {"path", "", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 2),
                    std::next(ps.begin(), 4),
                    "");
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "path/to/the/file.txt", "path/to/", {"path", "to", ""});
            check(f, "/path/to/the/file.txt", "/path/to/", {"path", "to", ""});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/the/file.txt", "test/the/file.txt", {"test", "the", "file.txt"});
            check(f, "/path/to/the/file.txt", "/test/the/file.txt", {"test", "the", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/the/file.txt", "path/test/file.txt", {"path", "test", "file.txt"});
            check(f, "/path/to/the/file.txt", "/path/test/file.txt", {"path", "test", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 2),
                    std::next(ps.begin(), 4),
                    "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/the/file.txt", "path/to/test", {"path", "to", "test"});
            check(f, "/path/to/the/file.txt", "/path/to/test", {"path", "to", "test"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    BIGSTR);
                BOOST_TEST_EQ(*it, BIGSTR);
            };
            check(f, "path/to/the/file.txt", "path/" BIGSTR "/file.txt",
                {"path", BIGSTR, "file.txt"});
        }

        //
        // replace(iterator, iterator. initializer_list)
        // replace(iterator, iterator. FwdIt, FwdIt)
        //

        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    { "t", "u", "v" });
                BOOST_TEST_EQ(*it, "t");
            };
            auto const g = [](segments_ref ps)
            {
                auto it = replace(ps,
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    { "t", "u", "v" });
                ignore_unused(it);
            };
            check(f, g, "path/to/the/file.txt", "t/u/v/the/file.txt", {"t", "u", "v", "the", "file.txt"});
            check(f, g, "/path/to/the/file.txt", "/t/u/v/the/file.txt", {"t", "u", "v", "the", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    { "t", "u", "v" });
                BOOST_TEST_EQ(*it, "t");
            };
            auto const g = [](segments_ref ps)
            {
                auto it = replace(ps,
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    { "t", "u", "v" });
                ignore_unused(it);
            };
            check(f, g, "path/to/the/file.txt", "path/t/u/v/file.txt", {"path", "t", "u", "v", "file.txt"});
            check(f, g, "/path/to/the/file.txt", "/path/t/u/v/file.txt", {"path", "t", "u", "v", "file.txt"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 2),
                    std::next(ps.begin(), 4),
                    { "t", "u", "v" });
                BOOST_TEST_EQ(*it, "t");
            };
            auto const g = [](segments_ref ps)
            {
                auto it = replace(ps,
                    std::next(ps.begin(), 2),
                    std::next(ps.begin(), 4),
                    { "t", "u", "v" });
                ignore_unused(it);
            };
            check(f, g, "path/to/the/file.txt", "path/to/t/u/v", {"path", "to", "t", "u", "v"});
            check(f, g, "/path/to/the/file.txt", "/path/to/t/u/v", {"path", "to", "t", "u", "v"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    { BIGSTR, BIGSTR, BIGSTR });
                BOOST_TEST_EQ(*it, BIGSTR);
            };
            auto const g = [](segments_ref ps)
            {
                auto it = replace(ps,
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    { BIGSTR, BIGSTR, BIGSTR });
                ignore_unused(it);
            };
            check(f, g, "path/to/the/file.txt",
                "path/" BIGSTR "/" BIGSTR "/" BIGSTR "/file.txt",
                {"path", BIGSTR, BIGSTR, BIGSTR, "file.txt"});
        }

        //
        // push_back
        //

        {
            auto const f = [](segments_ref ps)
            {
                ps.push_back("");
            };
            check(f, "",    "./",   {""});
            check(f, "/",   "/./",  {""});
            check(f, "./",  ".//",  {"", ""});
            check(f, "/./", "/.//", {"", ""});
        }
        {
            auto const f = [](segments_ref ps)
            {
                ps.push_back("/");
            };
            check(f, "",  "%2F",  {"/"});
            check(f, "/", "/%2F", {"/"});
        }
        {
            auto const f = [](segments_ref ps)
            {
                ps.push_back(":");
            };
            check(f, "",  "%3A",  {":"});
            check(f, "/", "/:", {":"});
        }

        //
        // pop_back
        //
        {
            auto const f = [](segments_ref ps)
            {
                ps.pop_back();
            };
            check(f, "/path/to/file.txt",  "/path/to",  {"path", "to"});
            check(f, "/path/to/",  "/path/to",  {"path", "to"});
            check(f, ".//",  "./",  {""});
            check(f, "/.//",  "/./",  {""});
            check(f, "x://y//",  "/",  {""});
            check(f, "x://y/.//",  "/./",  {""});
            check(f, "x://y/.///",  "/.//",  {"", ""});
        }
    }

    static
    void
    testJavadocs()
    {
        // {class}
        {
    url u( "/path/to/file.txt" );

    segments_ref ps = u.segments();

    ignore_unused(ps);
        }

        // operator=(initializer_list)
        {
        url u;

        u.segments() = { "path", "to", "file.txt" };
        }

        // url()
        {
        url u( "?key=value" );

        assert( &u.segments().url() == &u );

        ignore_unused(u);
        }

        // assign(initializer_list)
        {
        url u;
         
        u.segments().assign( { "path", "to", "file.txt" } );
        }

        // insert(iterator, initializer_list)
        {
        url u( "/file.txt" );

        u.segments().insert( u.segments().begin(), { "path", "to" } );
        }
    }

    //--------------------------------------------

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
    segments_ref_test,
    "boost.url.segments_ref");

} // urls
} // boost
