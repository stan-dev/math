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
#include <boost/url/segments_encoded_ref.hpp>

#include <boost/url/parse.hpp>
#include <boost/url/url.hpp>
#include <boost/static_assert.hpp>
#include <boost/core/ignore_unused.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

#ifdef assert
#undef assert
#endif
#define assert BOOST_TEST

BOOST_STATIC_ASSERT(
    ! std::is_default_constructible<
        segments_encoded_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_constructible<
        segments_encoded_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_copy_assignable<
        segments_encoded_ref>::value);

BOOST_STATIC_ASSERT(
    std::is_default_constructible<
        segments_encoded_ref::iterator>::value);

struct segments_encoded_ref_test
{
    // check that modification produces
    // the string and correct sequence
    static
    void
    check(
        void(*f)(segments_encoded_ref),
        string_view s0,
        string_view s1,
        std::initializer_list<
            string_view> init)
    {
        auto rv = parse_uri_reference(s0);
        if(! BOOST_TEST(rv.has_value()))
            return;
        url u = *rv;
        segments_encoded_ref ps(u.encoded_segments());
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
        void(*f1)(segments_encoded_ref),
        void(*f2)(segments_encoded_ref),
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
        segments_encoded_ref& ps,
        std::initializer_list<
            string_view> init)
    {
        ps.assign(init.begin(), init.end());
    }

    static
    auto
    insert(
        segments_encoded_ref& ps,
        segments_encoded_ref::iterator before,
        std::initializer_list<
            string_view> init) ->
        segments_encoded_ref::iterator
    {
        return ps.insert(before,
            init.begin(), init.end());
    }

    static
    auto
    replace(
        segments_encoded_ref& ps,
        segments_encoded_ref::iterator from,
        segments_encoded_ref::iterator to,
        std::initializer_list<
            string_view> init) ->
        segments_encoded_ref::iterator
    {
        return ps.replace(from, to,
            init.begin(), init.end());
    }

    //--------------------------------------------

    void
    testSpecial()
    {
        // segments_encoded_ref(segments_encoded_ref)
        {
            url u("/index.htm");
            segments_encoded_ref ps0 = u.encoded_segments();
            segments_encoded_ref ps1(ps0);
            BOOST_TEST_EQ(&ps0.url(), &ps1.url());
            BOOST_TEST_EQ(
                ps0.url().buffer().data(),
                ps1.url().buffer().data());
        }

        // operator=(segments_encoded_ref)
        {
            url u1("/index.htm");
            url u2("/path/to/file.txt");
            segments_encoded_ref ps1 = u1.encoded_segments();
            segments_encoded_ref ps2 = u2.encoded_segments();
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

        // operator=(segments_encoded_view)
        {
            url u1("/index.htm");
            url_view u2("/path/to/file.txt");
            segments_encoded_ref ps1 =
                u1.encoded_segments();
            segments_encoded_view ps2 =
                u2.encoded_segments();
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
            u.encoded_segments() = { "path", "to%3F", "file#" };
            BOOST_TEST_EQ(
                u.encoded_path(), "path/to%3F/file%23");
        }

        // operator segments_encoded_view()
        {
            url u;
            u.encoded_segments() = { "path", "to%3F", "file#" };
            segments_encoded_view ps = u.encoded_segments();
            auto it = ps.begin();
            BOOST_TEST_EQ(*it++, "path");
            BOOST_TEST_EQ(*it++, "to%3F");
            BOOST_TEST_EQ(*it++, "file%23");
        }

        // operator segments_encoded_view()
        {
            url u;
            u.encoded_segments() = { "x:y", "a:b" };
            segments_encoded_view ps = u.encoded_segments();
            auto it = ps.begin();
            BOOST_TEST_CSTR_EQ(*it++, "x%3Ay");
            BOOST_TEST_CSTR_EQ(*it++, "a:b");
        }
    }

    void
    testObservers()
    {
        // url()
        {
            url u0( "/" );
            url u1( "/" );
            BOOST_TEST_EQ(
                &u0.encoded_segments().url(), &u0);
            BOOST_TEST_EQ(
                &u1.encoded_segments().url(), &u1);
            BOOST_TEST_NE(
                &u0.encoded_segments().url(),
                &u1.encoded_segments().url());
        }
    }

    void
    testModifiers()
    {
        //
        // clear()
        //

        {
            auto const f = [](segments_encoded_ref ps)
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
            auto const f = [](segments_encoded_ref ps)
            {
                ps.assign({ "path", "to%23", "file.txt?" });

                // invalid percent-escape
                BOOST_TEST_THROWS(ps.assign({ "%" }), system_error);
            };
            auto const g = [](segments_encoded_ref ps)
            {
                assign(ps, { "path", "to%23", "file.txt?" });

                // invalid percent-escape
                BOOST_TEST_THROWS(assign(ps, { "%" }), system_error);
            };
            check(f, g, "", "path/to%23/file.txt%3F", {"path", "to%23", "file.txt%3F"});
            check(f, g, "/", "/path/to%23/file.txt%3F", {"path", "to%23", "file.txt%3F"});
            check(f, g, "/index.htm", "/path/to%23/file.txt%3F", {"path", "to%23", "file.txt%3F"});
            check(f, g, "index.htm", "path/to%23/file.txt%3F", {"path", "to%23", "file.txt%3F"});
            check(f, g, "/path/to/file.txt", "/path/to%23/file.txt%3F", {"path", "to%23", "file.txt%3F"});
            check(f, g, "Program%20Files", "path/to%23/file.txt%3F", {"path", "to%23", "file.txt%3F"});
        }

        //
        // insert(iterator, pct_string_view)
        //

        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.insert(ps.begin(), "");
                BOOST_TEST_EQ(*it, "");

                // invalid percent-escape
                BOOST_TEST_THROWS(ps.insert(ps.begin(),
                    "%%"), system_error);
            };
            check(f, "", "./", {""});
            check(f, "/", "/./", {""});
            check(f, "/index.htm", "/.//index.htm", {"", "index.htm"});
            check(f, "index.htm", ".//index.htm", {"", "index.htm"});
            check(f, "path/to/file.txt", ".//path/to/file.txt", {"", "path", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/.//path/to/file.txt", {"", "path", "to", "file.txt"});
            check(f, "Program%20Files", ".//Program%20Files", {"", "Program%20Files"});
            check(f, "x:", "./", {""});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.insert(ps.begin(), "my seg%23");
                BOOST_TEST_EQ(*it, "my%20seg%23");
            };
            check(f, "", "my%20seg%23", {"my%20seg%23"});
            check(f, "/", "/my%20seg%23", {"my%20seg%23"});
            check(f, "/index.htm", "/my%20seg%23/index.htm", {"my%20seg%23", "index.htm"});
            check(f, "index.htm", "my%20seg%23/index.htm", {"my%20seg%23", "index.htm"});
            check(f, "path/to/file.txt", "my%20seg%23/path/to/file.txt", {"my%20seg%23", "path", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/my%20seg%23/path/to/file.txt", {"my%20seg%23", "path", "to", "file.txt"});
            check(f, "Program%20Files", "my%20seg%23/Program%20Files", {"my%20seg%23", "Program%20Files"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.insert(std::next(ps.begin(), 1), "my%20seg?");
                BOOST_TEST_EQ(*it, "my%20seg%3F");
            };
            check(f, "path/to/file.txt", "path/my%20seg%3F/to/file.txt", {"path", "my%20seg%3F", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/path/my%20seg%3F/to/file.txt", {"path", "my%20seg%3F", "to", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.insert(ps.end(), "my%20seg[");
                BOOST_TEST_EQ(*it, "my%20seg%5B");
            };
            check(f, "", "my%20seg%5B", {"my%20seg%5B"});
            check(f, "/", "/my%20seg%5B", {"my%20seg%5B"});
            check(f, "/index.htm", "/index.htm/my%20seg%5B", {"index.htm", "my%20seg%5B"});
            check(f, "index.htm", "index.htm/my%20seg%5B", {"index.htm", "my%20seg%5B"});
            check(f, "path/to/file.txt", "path/to/file.txt/my%20seg%5B", {"path", "to", "file.txt", "my%20seg%5B"});
            check(f, "/path/to/file.txt", "/path/to/file.txt/my%20seg%5B", {"path", "to", "file.txt", "my%20seg%5B"});
            check(f, "Program%20Files", "Program%20Files/my%20seg%5B", {"Program%20Files", "my%20seg%5B"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
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

        //
        // insert(iterator, initializer_list)
        // insert(iterator, FwdIt, FwdIt)
        //

        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.insert(ps.begin(), { "u#", "v%20" });
                BOOST_TEST_EQ(*it, "u%23");

                // invalid percent-escape
                BOOST_TEST_THROWS(ps.insert(ps.begin(),
                    { "%F" }), system_error);
            };
            auto const g = [](segments_encoded_ref ps)
            {
                auto it = insert(ps, ps.begin(), { "u#", "v%20" });
                BOOST_TEST_EQ(*it, "u%23");

                // invalid percent-escape
                BOOST_TEST_THROWS(insert(ps, ps.begin(),
                    { "%F" }), system_error);
            };
            check(f, g, "", "u%23/v%20", {"u%23", "v%20"});
            check(f, g, "/", "/u%23/v%20", {"u%23", "v%20"});
            check(f, g, "/index.htm", "/u%23/v%20/index.htm", {"u%23", "v%20", "index.htm"});
            check(f, g, "index.htm", "u%23/v%20/index.htm", {"u%23", "v%20", "index.htm"});
            check(f, g, "path/to/file.txt", "u%23/v%20/path/to/file.txt", {"u%23", "v%20", "path", "to", "file.txt"});
            check(f, g, "/path/to/file.txt", "/u%23/v%20/path/to/file.txt", {"u%23", "v%20", "path", "to", "file.txt"});
            check(f, g, "Program%20Files", "u%23/v%20/Program%20Files", {"u%23", "v%20", "Program%20Files"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.insert(ps.begin(), { "", "" });
                BOOST_TEST_EQ(*it, "");
            };
            auto const g = [](segments_encoded_ref ps)
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

        //
        // erase(iterator)
        //

        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.erase(std::next(ps.begin(), 0));
                BOOST_TEST_EQ(*it, ps.front());
            };
            check(f, "path/to/file.txt", "to/file.txt", {"to", "file.txt"});
            check(f, "/path/to/file.txt", "/to/file.txt", {"to", "file.txt"});
            check(f, "//x/y/", "/./", {""});
            check(f, "/x/", "/./", {""});
            check(f, "x/", "./", {""});
            check(f, "x:.//", "./", {""});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.erase(std::next(ps.begin(), 1));
                BOOST_TEST_EQ(*it, "file.txt");
            };
            check(f, "path/to/file.txt", "path/file.txt", {"path", "file.txt"});
            check(f, "/path/to/file.txt", "/path/file.txt", {"path", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.erase(std::next(ps.begin(), 2));
                BOOST_TEST_EQ(it, ps.end());
            };
            check(f, "path/to/file.txt", "path/to", {"path", "to"});
            check(f, "/path/to/file.txt", "/path/to", {"path", "to"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.erase(std::next(ps.begin(), 1));
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "x://y///", "//", {"", ""});
            check(f, ".///", ".//", {"", ""});
        }

        //
        // erase(iterator, iterator)
        //

        {
            auto const f = [](segments_encoded_ref ps)
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
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.erase(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3));
                BOOST_TEST_EQ(*it, ps.back());
            };
            check(f, "path/to/the/file.txt", "path/file.txt", {"path", "file.txt"});
            check(f, "/path/to/the/file.txt", "/path/file.txt", {"path", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
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
        // replace(iterator, pct_string_view)
        //

        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 0), "");
                BOOST_TEST_EQ(*it, "");

                // invalid percent escape
                BOOST_TEST_THROWS(ps.replace(std::next(ps.begin(), 0),
                    "00%"), system_error);
            };
            check(f, "path/to/file.txt", ".//to/file.txt", {"", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/.//to/file.txt", {"", "to", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 1), "");
                BOOST_TEST_EQ(*it, "");
            };
            check(f, "path/to/file.txt", "path//file.txt", {"path", "", "file.txt"});
            check(f, "/path/to/file.txt", "/path//file.txt", {"path", "", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 0), "te%20[");
                BOOST_TEST_EQ(*it, "te%20%5B");
            };
            check(f, "path/to/file.txt", "te%20%5B/to/file.txt", {"te%20%5B", "to", "file.txt"});
            check(f, "/path/to/file.txt", "/te%20%5B/to/file.txt", {"te%20%5B", "to", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 1), "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/file.txt", "path/test/file.txt", {"path", "test", "file.txt"});
            check(f, "/path/to/file.txt", "/path/test/file.txt", {"path", "test", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(std::next(ps.begin(), 2), "test");
                BOOST_TEST_EQ(*it, "test");
            };
            check(f, "path/to/file.txt", "path/to/test", {"path", "to", "test"});
            check(f, "/path/to/file.txt", "/path/to/test", {"path", "to", "test"});
        }

        //
        // replace(iterator, iterator, pct_string_view)
        //

        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    "");
                BOOST_TEST_EQ(*it, "");

                // invalid percent escape
                BOOST_TEST_THROWS(ps.replace(
                    std::next(ps.begin(), 0), std::next(ps.begin(), 2),
                    "0%"), system_error);
            };
            check(f, "path/to/the/file.txt", ".//the/file.txt", {"", "the", "file.txt"});
            check(f, "/path/to/the/file.txt", "/.//the/file.txt", {"", "the", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
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
            auto const f = [](segments_encoded_ref ps)
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
            auto const f = [](segments_encoded_ref ps)
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
            auto const f = [](segments_encoded_ref ps)
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
            auto const f = [](segments_encoded_ref ps)
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

        //
        // replace(iterator, iterator. initializer_list)
        // replace(iterator, iterator. FwdIt, FwdIt)
        //

        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    { "t", "u %3F", "v" });
                BOOST_TEST_EQ(*it, "t");

                // invalid percent escape
                BOOST_TEST_THROWS(ps.replace(
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    { "x", "%FG" }), system_error);
            };
            auto const g = [](segments_encoded_ref ps)
            {
                auto it = replace(ps,
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    { "t", "u %3F", "v" });
                BOOST_TEST_EQ(*it, "t");

                // invalid percent escape
                BOOST_TEST_THROWS(replace(ps,
                    std::next(ps.begin(), 0),
                    std::next(ps.begin(), 2),
                    { "x", "%" }), system_error);
            };
            check(f, g, "path/to/the/file.txt", "t/u%20%3F/v/the/file.txt", {"t", "u%20%3F", "v", "the", "file.txt"});
            check(f, g, "/path/to/the/file.txt", "/t/u%20%3F/v/the/file.txt", {"t", "u%20%3F", "v", "the", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    { "t", "u", "v" });
                BOOST_TEST_EQ(*it, "t");
            };
            auto const g = [](segments_encoded_ref ps)
            {
                auto it = replace(ps,
                    std::next(ps.begin(), 1),
                    std::next(ps.begin(), 3),
                    { "t", "u", "v" });
                BOOST_TEST_EQ(*it, "t");
            };
            check(f, g, "path/to/the/file.txt", "path/t/u/v/file.txt", {"path", "t", "u", "v", "file.txt"});
            check(f, g, "/path/to/the/file.txt", "/path/t/u/v/file.txt", {"path", "t", "u", "v", "file.txt"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                auto it = ps.replace(
                    std::next(ps.begin(), 2),
                    std::next(ps.begin(), 4),
                    { "t", "u", "v" });
                BOOST_TEST_EQ(*it, "t");
            };
            auto const g = [](segments_encoded_ref ps)
            {
                auto it = replace(ps,
                    std::next(ps.begin(), 2),
                    std::next(ps.begin(), 4),
                    { "t", "u", "v" });
                BOOST_TEST_EQ(*it, "t");
            };
            check(f, g, "path/to/the/file.txt", "path/to/t/u/v", {"path", "to", "t", "u", "v"});
            check(f, g, "/path/to/the/file.txt", "/path/to/t/u/v", {"path", "to", "t", "u", "v"});
        }

        //
        // push_back
        //

        {
            auto const f = [](segments_encoded_ref ps)
            {
                ps.push_back("");

                // invalid percent-escape
                BOOST_TEST_THROWS(ps.push_back("%"), system_error);
            };
            check(f, "",    "./",   {""});
            check(f, "/",   "/./",  {""});
            check(f, "./",  ".//",  {"", ""});
            check(f, "/./", "/.//", {"", ""});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                ps.push_back("/");
            };
            check(f, "",  "%2F",  {"%2F"});
            check(f, "/", "/%2F", {"%2F"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                ps.push_back(":");
            };
            check(f, "",  "%3A",  {"%3A"});
            check(f, "/", "/:", {":"});
        }

        //
        // pop_back
        //
        {
            auto const f = [](segments_encoded_ref ps)
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

    void
    testEditSegments()
    {
    /*  Legend

        '#' 0x23    '/' 0x2f
        '%' 0x25    ':' 0x3a
        '.' 0x2e    '?' 0x3f
    */
        {
            auto const f = [](segments_encoded_ref ps)
            {
                ps.push_back("");
            };
            check(f, "",    "./",   {""});
            check(f, "/",   "/./",  {""});
            check(f, "./",  ".//",  {"", ""});
            check(f, "/./", "/.//", {"", ""});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                ps.push_back("/");
            };
            check(f, "",  "%2F",  {"%2F"});
            check(f, "/", "/%2F", {"%2F"});
        }
        {
            auto const f = [](segments_encoded_ref ps)
            {
                ps.push_back(":");
            };
            check(f, "",  "%3A",  {"%3A"});
            check(f, "/", "/:", {":"});
        }
    }

    void
    testJavadocs()
    {
        // {class}
        {
    url u( "/path/to/file.txt" );

    segments_encoded_ref ps = u.encoded_segments();

    ignore_unused(ps);
        }

        // operator=(initializer_list)
        {
        url u;

        u.encoded_segments() = {"path", "to", "file.txt"};
        }

        // url()
        {
        url u( "?key=value" );

        assert( &u.encoded_segments().url() == &u );
        }
    }

    void
    run()
    {
        testSpecial();
        testObservers();
        testModifiers();
        testEditSegments();
        testJavadocs();
    }
};

TEST_SUITE(
    segments_encoded_ref_test,
    "boost.url.segments_encoded_ref");

} // urls
} // boost
