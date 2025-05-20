//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Antony Polukhin, 2014-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <boost/lexical_cast/try_lexical_convert.hpp>

#include <boost/core/lightweight_test.hpp>

using namespace boost::conversion;

void try_uncommon_cases()
{
    std::string sres;
    const bool res1 = try_lexical_convert(std::string("Test string"), sres);
    BOOST_TEST(res1);
    BOOST_TEST_EQ(sres, "Test string");

    volatile int vires;
    const bool res2 = try_lexical_convert(100, vires);
    BOOST_TEST(res2);
    BOOST_TEST_EQ(vires, 100);

    const bool res3 = try_lexical_convert("Test string", sres);
    BOOST_TEST(res3);
    BOOST_TEST_EQ(sres, "Test string");

    const bool res4 = try_lexical_convert("Test string", sizeof("Test string") - 1, sres);
    BOOST_TEST(res4);
    BOOST_TEST_EQ(sres, "Test string");

    int ires;
    BOOST_TEST(!try_lexical_convert("Test string", ires));
    BOOST_TEST(!try_lexical_convert(1.1, ires));
    BOOST_TEST(!try_lexical_convert(-1.9, ires));
    BOOST_TEST(!try_lexical_convert("1.1", ires));
    BOOST_TEST(!try_lexical_convert("1000000000000000000000000000000000000000", ires));
}


void try_common_cases()
{
    int ires = 0;
    const bool res1 = try_lexical_convert(std::string("100"), ires);
    BOOST_TEST(res1);
    BOOST_TEST_EQ(ires, 100);

    ires = 0;
    const bool res2 = try_lexical_convert("-100", ires);
    BOOST_TEST(res2);
    BOOST_TEST_EQ(ires, -100);

    float fres = 1.0f;
    const bool res3 = try_lexical_convert("0.0", fres);
    BOOST_TEST(res3);
    BOOST_TEST_EQ(fres, 0.0f);

    fres = 1.0f;
    const bool res4 = try_lexical_convert("0.0", sizeof("0.0") - 1, fres);
    BOOST_TEST(res4);
    BOOST_TEST_EQ(fres, 0.0f);
}

int main()
{
    try_uncommon_cases();
    try_common_cases();

    return boost::report_errors();
}
