//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#ifndef BOOST_URL_SOURCE
#define BOOST_URL_SOURCE
#endif

#include "router.hpp"

#include "test_suite.hpp"

namespace boost {
namespace urls {

struct router_test
{
    static
    void
    testPatterns()
    {
        // literal
        good("user", "user");
        good("user/view", "user/view");

        // empty
        good("", "/");
        good("/", "/");
        good("", "");
        good("/", "");

        // replace segment match
        good("user/{}", "user/johndoe", {"johndoe"});
        good("user/{name}", "user/johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good("user/{id}", "user/123", {"123"}, {{"id", "123"}});
        good("user/{name}/{op}", "user/johndoe/r", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        good("user/{name}/op/{op}", "user/johndoe/op/r", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});

        // path dot segments
        good("user/{name}", "user/././johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good("user/{name}", "user/b/../johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good("user/b", "user/c/../b");
        good("user/b", "../a/user/c/../b");

        // template dot segments
        good("user/././{name}", "user/johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good("user/b/../{name}", "user/johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good("user/c/../b", "user/b");
        good("user/{name}/../{name}", "user/johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good("../a/user/c/../b", "user/b");

        // pct-encoded segments
        good("%75ser/{name}", "user/johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good("user/{name}", "%75ser/johndoe", {"johndoe"}, {{"name", "johndoe"}});

        // empty segments
        good("user//{name}", "user//johndoe", {"johndoe"}, {{"name", "johndoe"}});
        bad("user//{name}", "user/johndoe");

        // optional segment match
        good("user/{name}/{op?}", "user/johndoe", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op?}", "user/johndoe/", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op?}", "user/johndoe/r", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        bad("user/{name}/{op?}", "user/johndoe/r/r");
        good("user/{name}/{op?}/b", "user/johndoe/r/b", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        good("user/{name}/{op?}/b", "user/johndoe/b", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op?}/{op2?}", "user/johndoe", {"johndoe", "", ""}, {{"name", "johndoe"}, {"op", ""}, {"op2", ""}});
        good("user/{name}/{op?}/{op2?}", "user/johndoe/", {"johndoe", "", ""}, {{"name", "johndoe"}, {"op", ""}, {"op2", ""}});
        good("user/{name}/{op?}/{op2?}", "user/johndoe/r", {"johndoe", "r", ""}, {{"name", "johndoe"}, {"op", "r"}, {"op2", ""}});
        good("user/{name}/{op?}/{op2?}/b", "user/johndoe/b", {"johndoe", "", ""}, {{"name", "johndoe"}, {"op", ""}, {"op2", ""}});
        good("user/{name}/{op?}/{op2?}/b", "user/johndoe//b", {"johndoe", "", ""}, {{"name", "johndoe"}, {"op", ""}, {"op2", ""}});
        good("user/{name}/{op?}/{op2?}/b", "user/johndoe/r/b", {"johndoe", "r", ""}, {{"name", "johndoe"}, {"op", "r"}, {"op2", ""}});

        // optional segment match + dot segments
        good("user/{name}/{op?}", "user/johndoe/.", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op?}", "user/johndoe/./r", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        bad("user/{name}/{op?}", "user/johndoe/..");
        good("user/{name}/{op?}", "user/johndoe/../maria", {"maria", ""}, {{"name", "maria"}, {"op", ""}});
        good("user/{name}/{op?}", "user/johndoe/../maria/./r", {"maria", "r"}, {{"name", "maria"}, {"op", "r"}});
        good("user/{name}/{op?}", "user/johndoe/r/../s", {"johndoe", "s"}, {{"name", "johndoe"}, {"op", "s"}});

        // star segment
        good("user/{name}/{op*}", "user/johndoe", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op*}", "user/johndoe/", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op*}", "user/johndoe/r", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        good("user/{name}/{op*}", "user/johndoe/r/r1", {"johndoe", "r/r1"}, {{"name", "johndoe"}, {"op", "r/r1"}});
        good("user/{name}/{op*}", "user/johndoe/r/r1/r2", {"johndoe", "r/r1/r2"}, {{"name", "johndoe"}, {"op", "r/r1/r2"}});
        good("user/{name}/{op*}/{op2*}", "user/johndoe", {"johndoe", "", ""}, {{"name", "johndoe"}, {"op", ""}, {"op2", ""}});
        good("user/{name}/{op*}/{op2*}", "user/johndoe/", {"johndoe", "", ""}, {{"name", "johndoe"}, {"op", ""}, {"op2", ""}});
        good("user/{name}/{op*}/{op2*}", "user/johndoe/r", {"johndoe", "r", ""}, {{"name", "johndoe"}, {"op", "r"}, {"op2", ""}});
        good("user/{name}/{op*}/b", "user/johndoe/b", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op*}/b", "user/johndoe/r/b", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        good("user/{name}/{op*}/b", "user/johndoe/r/r1/b", {"johndoe", "r/r1"}, {{"name", "johndoe"}, {"op", "r/r1"}});
        good("user/{name}/{op*}/b", "user/johndoe/r/r1/r2/b", {"johndoe", "r/r1/r2"}, {{"name", "johndoe"}, {"op", "r/r1/r2"}});
        good("user/{name}/{op*}/b", "user/johndoe/b", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});

        // star segment match + dot segments
        good("user/{name}/{op*}", "user/johndoe/../r/../s", {"s", ""}, {{"name", "s"}, {"op", ""}});
        bad("user/{name}/{op*}", "user/johndoe/../../r/s");
        good("user/{name}/{op*}", "user/johndoe/../../user/s", {"s", ""}, {{"name", "s"}, {"op", ""}});
        good("user/{name}/{op*}", "user/johndoe/r/../../../user/s", {"s", ""}, {{"name", "s"}, {"op", ""}});
        good("user/{name}/{op*}", "user/johndoe/r/../s", {"johndoe", "r/../s"}, {{"name", "johndoe"}, {"op", "r/../s"}});
        good("user/{name}/{op*}", "user/johndoe/r/../s/..", {"johndoe", "r/../s/.."}, {{"name", "johndoe"}, {"op", "r/../s/.."}});
        bad("user/{name}/{op*}", "user/johndoe/r/../s/../..");

        // plus segment
        bad("user/{name}/{op+}", "user/johndoe");
        good("user/{name}/{op+}", "user/johndoe/", {"johndoe", ""}, {{"name", "johndoe"}, {"op", ""}});
        good("user/{name}/{op+}", "user/johndoe/r", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        good("user/{name}/{op+}", "user/johndoe/r/r1", {"johndoe", "r/r1"}, {{"name", "johndoe"}, {"op", "r/r1"}});
        good("user/{name}/{op+}", "user/johndoe/r/r1/r2", {"johndoe", "r/r1/r2"}, {{"name", "johndoe"}, {"op", "r/r1/r2"}});
        bad("user/{name}/{op+}/{op2+}", "user/johndoe");
        bad("user/{name}/{op+}/{op2+}", "user/johndoe/");
        bad("user/{name}/{op+}/{op2+}", "user/johndoe/r");
        bad("user/{name}/{op+}/b", "user/johndoe/b");
        good("user/{name}/{op+}/b", "user/johndoe/r/b", {"johndoe", "r"}, {{"name", "johndoe"}, {"op", "r"}});
        good("user/{name}/{op+}/b", "user/johndoe/r/r1/b", {"johndoe", "r/r1"}, {{"name", "johndoe"}, {"op", "r/r1"}});
        good("user/{name}/{op+}/b", "user/johndoe/r/r1/r2/b", {"johndoe", "r/r1/r2"}, {{"name", "johndoe"}, {"op", "r/r1/r2"}});
        bad("user/{name}/{op+}/b", "user/johndoe/b");

        // plus segment match + dot segments
        bad("user/{name}/{op+}", "user/johndoe/../r/../s");
        good("user/{name}/{op+}", "user/johndoe/../r/../s/", {"s", ""}, {{"name", "s"}, {"op", ""}});
        bad("user/{name}/{op+}", "user/johndoe/../../r/s");
        bad("user/{name}/{op+}", "user/johndoe/../../user/s");
        good("user/{name}/{op+}", "user/johndoe/../../user/s/", {"s", ""}, {{"name", "s"}, {"op", ""}});
        good("user/{name}/{op+}", "user/johndoe/../../user/s/r", {"s", "r"}, {{"name", "s"}, {"op", "r"}});
        bad("user/{name}/{op+}", "user/johndoe/r/../../../user/s");
        good("user/{name}/{op+}", "user/johndoe/r/../../../user/s/", {"s", ""}, {{"name", "s"}, {"op", ""}});
        good("user/{name}/{op+}", "user/johndoe/r/../s", {"johndoe", "s"}, {{"name", "johndoe"}, {"op", "s"}});
        bad("user/{name}/{op+}", "user/johndoe/r/../s/..");
        bad("user/{name}/{op+}", "user/johndoe/r/../s/../..");

        // segment tree
        good({"user/view", "user/{}"}, 1, "user/johndoe", {"johndoe"});
        good({"user/view", "user/{}"}, 0, "user/view", {});
        good({"user/{}", "user/view"}, 1, "user/view", {});
        good({"user/view", "user/{file*}"}, 1, "user/path/to/file.txt", {"path/to/file.txt"}, {{"file", "path/to/file.txt"}});
        good({"user/view", "user/{name}", "user/{file+}"}, 1, "user/johndoe", {"johndoe"}, {{"name", "johndoe"}});
        good({"user/view", "user/{name}", "user/{file+}"}, 2, "user/path/to/file.txt", {"path/to/file.txt"}, {{"file", "path/to/file.txt"}});
        good({"a/a", "a/b", "a/c", "a/d", "a/e", "a/f", "a/g", "a/h", "a/i", "a/j", "a/k"}, 0, "a/a", {}, {});

        // above root
        bad("user/../..");

        // neither segment nor field
        bad("user/{", "user/johndoe");
    }

    static
    void
    good(
        std::initializer_list<core::string_view> patterns,
        std::size_t match_idx,
        core::string_view request,
        std::initializer_list<core::string_view> ms = {},
        std::initializer_list<std::pair<core::string_view, core::string_view>> args = {})
    {
        router<int> r;
        int i = 0;
        for (auto p: patterns)
            r.insert(p, i++);
        matches m;
        int const* v;
        v = r.find(request, m);
        if (!BOOST_TEST(v))
            return;
        BOOST_TEST_EQ(*v, match_idx);

        BOOST_TEST_EQ(ms.size(), m.size());
        BOOST_TEST_EQ(ms.size() == 0, m.empty());

        {
            auto it0 = ms.begin();
            auto it1 = m.begin();
            auto end0 = ms.end();
            auto end1 = m.end();
            std::size_t j = 0;
            while (
                it0 != end0 &&
                it1 != end1)
            {
                BOOST_TEST_EQ(*it0, *it1++);
                BOOST_TEST_EQ(*it0, m[j]);
                BOOST_TEST_EQ(*it0++, m.at(j++));
            }
            BOOST_TEST_THROWS(m.at(j), std::out_of_range);
        }

        {
            auto it = args.begin();
            auto end = args.end();
            while (it != end)
            {
                BOOST_TEST(m.find(it->first) != m.end());
                BOOST_TEST_EQ(m[it->first], it->second);
                BOOST_TEST_EQ(m.at(it->first), it->second);
                ++it;
            }
            BOOST_TEST(m.find("XXX") == m.end());
            BOOST_TEST_THROWS(m.at("XXX"), std::out_of_range);
        }
    };

    static
    void
    good(
        core::string_view pattern,
        core::string_view request,
        std::initializer_list<core::string_view> matches = {},
        std::initializer_list<std::pair<core::string_view, core::string_view>> args = {})
    {
        return good({pattern}, 0, request, matches, args);
    };

    static
    void
    bad(
        core::string_view pattern,
        core::string_view request)
    {
        router<int> r;
        try
        {
            r.insert(pattern, 0);
        }
        catch (...)
        {
            BOOST_TEST(true);
            return;
        }
        matches m;
        int const* v;
        v = r.find(request, m);
        BOOST_TEST_NOT(v);
    };

    static
    void
    bad(core::string_view pattern)
    {
        router<int> r;
        BOOST_TEST_THROWS(r.insert(pattern, 0), system::system_error);
    };

    void
    run()
    {
        testPatterns();
    }
};

TEST_SUITE(router_test, "boost.url.router");

} // urls
} // boost
