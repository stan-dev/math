//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Antony Polukhin, 2011-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <boost/lexical_cast.hpp>

#include <boost/core/lightweight_test.hpp>

void test_typedefed_wchar_t(unsigned short)  // wchar_t is a typedef for unsigned short
{
    BOOST_TEST(boost::lexical_cast<int>(L'A') == 65);
    BOOST_TEST(boost::lexical_cast<int>(L'B') == 66);

    BOOST_TEST(boost::lexical_cast<wchar_t>(L"65") == 65);
    BOOST_TEST(boost::lexical_cast<wchar_t>(L"66") == 66);
}

template <class T>
void test_typedefed_wchar_t(T)
{
    BOOST_TEST(1);
}

void test_typedefed_wchar_t_runtime()
{
    test_typedefed_wchar_t(L'0');
}

int main()
{
    test_typedefed_wchar_t_runtime();

    return boost::report_errors();
}
