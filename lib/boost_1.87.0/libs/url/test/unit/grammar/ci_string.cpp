//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/grammar/ci_string.hpp>

#include <boost/url/string_view.hpp>
#include <boost/container/map.hpp>
#include <boost/unordered_map.hpp>
#include "test_suite.hpp"

#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>

namespace boost {
namespace urls {
namespace grammar {

class ascii_test
{
public:
    void
    testToUpperLower()
    {
        core::string_view s0 =
            "0123456789"
            "abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        core::string_view s1 =
            "0123456789"
            "abcdefghijklmnopqrstuvwxyz"
            "abcdefghijklmnopqrstuvwxyz";

        auto it0 = s0.data();
        auto it1 = s1.data();
        auto const end = it0 + s0.size();
        BOOST_TEST_EQ(s0.size(), s1.size());
        while(it0 != end)
        {
            BOOST_TEST(
                to_lower(*it0) == *it1);
            BOOST_TEST(
                to_lower(*it0) ==
                to_lower(*it1));
            BOOST_TEST(
                to_upper(*it0) ==
                to_upper(*it1));
            ++it0;
            ++it1;
        }
    }

    void
    testDigest()
    {
        BOOST_TEST_EQ(ci_digest(""),    ci_digest(""));
        BOOST_TEST_EQ(ci_digest("x"),   ci_digest("x"));
        BOOST_TEST_EQ(ci_digest("x"),   ci_digest("X"));
        BOOST_TEST_EQ(ci_digest("abc"), ci_digest("aBc"));

        BOOST_TEST_NE(ci_digest("ABC"), ci_digest("XYZ"));
        BOOST_TEST_NE(ci_digest("abc"), ci_digest("abcd"));

        BOOST_TEST_EQ(
            ci_hash{}("abc"), ci_hash{}("ABC"));
        BOOST_TEST_NE(
            ci_hash{}("xyz"), ci_hash{}(""));
    }

    void
    testIsEqual()
    {
        BOOST_TEST(ci_is_equal("", ""));
        BOOST_TEST(ci_is_equal("a", "a"));
        BOOST_TEST(ci_is_equal("A", "A"));
        BOOST_TEST(ci_is_equal("aBC", "aBc"));
        BOOST_TEST(! ci_is_equal("aBC", "aBcd"));
        BOOST_TEST(! ci_is_equal("aBC", ""));

        BOOST_TEST(ci_equal{}("abc", "ABC"));
        BOOST_TEST(! ci_equal{}("xz", "abc"));
    }

    void
    testIsLess()
    {
        BOOST_TEST(ci_is_less("a", "aa"));
        BOOST_TEST(ci_is_less("aa", "ab"));
        BOOST_TEST(ci_is_less("Aa", "aB"));
        BOOST_TEST(! ci_is_less("BB", "b"));

        BOOST_TEST(ci_less{}("a", "aa"));
        BOOST_TEST(! ci_less{}("xy", "z"));
    }

    void
    testCompare()
    {
        BOOST_TEST_EQ(ci_compare("a", "b"), -1);
        BOOST_TEST_EQ(ci_compare("A", "b"), -1);
        BOOST_TEST_EQ(ci_compare("a", "B"), -1);
        BOOST_TEST_EQ(ci_compare("A", "B"), -1);

        BOOST_TEST_EQ(ci_compare("b", "a"), 1);
        BOOST_TEST_EQ(ci_compare("b", "A"), 1);
        BOOST_TEST_EQ(ci_compare("B", "a"), 1);
        BOOST_TEST_EQ(ci_compare("B", "A"), 1);

        BOOST_TEST_EQ(ci_compare("a", "a"), 0);
        BOOST_TEST_EQ(ci_compare("a", "A"), 0);
        BOOST_TEST_EQ(ci_compare("A", "a"), 0);
        BOOST_TEST_EQ(ci_compare("A", "A"), 0);

        BOOST_TEST_EQ(ci_compare("a", "ab"),  -1);
        BOOST_TEST_EQ(ci_compare("ab", "a"),   1);
        BOOST_TEST_EQ(ci_compare("aa", "ab"), -1);
        BOOST_TEST_EQ(ci_compare("ba", "bb"), -1);

        BOOST_TEST_EQ(ci_compare("A", "ab"),  -1);
        BOOST_TEST_EQ(ci_compare("ab", "A"),   1);
        BOOST_TEST_EQ(ci_compare("Aa", "aB"), -1);
        BOOST_TEST_EQ(ci_compare("bA", "BB"), -1);
    }

    void
    run()
    {
        // javadocs
        {
            assert( to_lower( 'A' ) == 'a' );
            assert( to_upper( 'a' ) == 'A' );
            assert( ci_compare( "boost", "Boost" ) == 0 );
            assert( ci_is_equal( "Boost", "boost" ) );
            assert( ! ci_is_less( "Boost", "boost" ) );
            {
                boost::unordered_map< std::string, std::string, ci_hash, ci_equal > m1;

                std::unordered_map  < std::string, std::string, ci_hash, ci_equal > m2; // (since C++20)

                m1.emplace("key", "value");
                BOOST_TEST(m1.count("KEY") == 1);
            }
            {
                boost::container::map< std::string, std::string, ci_less > m1;

                std::map  < std::string, std::string, ci_less > m2; // (since C++14)

                m1.emplace("key", "value");
                BOOST_TEST(m1.count("KEY") == 1);
            }
        }

        testToUpperLower();
        testDigest();
        testIsEqual();
        testIsLess();
        testCompare();
    }
};

TEST_SUITE(
    ascii_test,
    "boost.url.grammar.ci_string");

} // grammar
} // urls
} // boost
