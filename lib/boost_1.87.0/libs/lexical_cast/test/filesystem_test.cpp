//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Antony Polukhin, 2013-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).
//
// Test lexical_cast usage with long filesystem::path. Bug 7704.

#include <boost/lexical_cast.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/filesystem/path.hpp>

void test_filesystem()
{
    boost::filesystem::path p;
    std::string s1 = "aaaaaaaaaaaaaaaaaaaaaaa";
    p = boost::lexical_cast<boost::filesystem::path>(s1);
    BOOST_TEST(!p.empty());
    BOOST_TEST_EQ(p, s1);
    p.clear();

    const char ab[] = "aaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
    p = boost::lexical_cast<boost::filesystem::path>(ab);
    BOOST_TEST(!p.empty());
    BOOST_TEST_EQ(p, ab);

    // Tests for
    // https://github.com/boostorg/lexical_cast/issues/25

    const char quoted_path[] = "\"/home/my user\"";
    p = boost::lexical_cast<boost::filesystem::path>(quoted_path);
    BOOST_TEST(!p.empty());
    const char unquoted_path[] = "/home/my user";
    BOOST_TEST_EQ(p, boost::filesystem::path(unquoted_path));

    // Converting back to std::string gives the initial string
    BOOST_TEST_EQ(boost::lexical_cast<std::string>(p), quoted_path);

    try {
        // Without quotes the path will have only `/home/my` in it.
        // `user` remains in the stream, so an exception must be thrown.
        p = boost::lexical_cast<boost::filesystem::path>(unquoted_path);
        BOOST_TEST(false);
    } catch (const boost::bad_lexical_cast& ) {
        // Exception is expected
    }
}

int main()
{
    test_filesystem();

    return boost::report_errors();
}

