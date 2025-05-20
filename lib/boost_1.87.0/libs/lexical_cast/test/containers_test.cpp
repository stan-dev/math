//  Testing boost::lexical_cast with boost::container::string.
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
#include <boost/container/string.hpp>

using namespace boost;

void testing_boost_containers_basic_string()
{
    BOOST_TEST("100" == lexical_cast<boost::container::string>("100"));
    BOOST_TEST(L"100" == lexical_cast<boost::container::wstring>(L"100"));

    BOOST_TEST("100" == lexical_cast<boost::container::string>(100));
    boost::container::string str("1000");
    BOOST_TEST(1000 == lexical_cast<int>(str));
}

#if defined(BOOST_NO_STRINGSTREAM) || defined(BOOST_NO_STD_WSTRING)
#define BOOST_LCAST_NO_WCHAR_T
#endif

void testing_boost_containers_string_std_string()
{
    std::string std_str("std_str");
    boost::container::string boost_str("boost_str");
    BOOST_TEST(boost::lexical_cast<std::string>(boost_str) == "boost_str");
    BOOST_TEST(boost::lexical_cast<boost::container::string>(std_str) == "std_str");

#ifndef BOOST_LCAST_NO_WCHAR_T
    std::wstring std_wstr(L"std_wstr");
    boost::container::wstring boost_wstr(L"boost_wstr");

    BOOST_TEST(boost::lexical_cast<std::wstring>(boost_wstr) == L"boost_wstr");
    BOOST_TEST(boost::lexical_cast<boost::container::wstring>(std_wstr) == L"std_wstr");

#endif

}

void testing_boost_containers_string_widening()
{
    const char char_array[] = "Test string";

#ifndef BOOST_LCAST_NO_WCHAR_T
    const wchar_t wchar_array[] = L"Test string";
    BOOST_TEST(boost::lexical_cast<boost::container::wstring>(char_array) == wchar_array);
#endif

#if !defined(BOOST_NO_CXX11_CHAR16_T) && !defined(BOOST_NO_CXX11_UNICODE_LITERALS) && defined(BOOST_STL_SUPPORTS_NEW_UNICODE_LOCALES)
    const char16_t char16_array[] = u"Test string";
    BOOST_TEST(boost::lexical_cast<boost::container::basic_string<char16_t> >(char_array) == char16_array);
#endif

#if !defined(BOOST_NO_CXX11_CHAR32_T) && !defined(BOOST_NO_CXX11_UNICODE_LITERALS) && defined(BOOST_STL_SUPPORTS_NEW_UNICODE_LOCALES)
    const char32_t char32_array[] = U"Test string";
    BOOST_TEST(boost::lexical_cast<boost::container::basic_string<char32_t> >(char_array) == char32_array);
#endif
}

int main()
{
    testing_boost_containers_basic_string();
    testing_boost_containers_string_std_string();
    testing_boost_containers_string_widening();

    return boost::report_errors();
}
