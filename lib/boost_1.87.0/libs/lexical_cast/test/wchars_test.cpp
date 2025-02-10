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

#if defined(BOOST_NO_STRINGSTREAM) || defined(BOOST_NO_STD_WSTRING)
#define BOOST_LCAST_NO_WCHAR_T
#endif

template <class CharT>
void test_impl(const CharT* wc_arr)
{
    typedef CharT                       wide_char;
    typedef std::basic_string<CharT>    wide_string;
    const char c_arr[]            = "Test array of chars";
    const unsigned char uc_arr[]  = "Test array of chars";
    const signed char sc_arr[]    = "Test array of chars";

    // Following tests depend on realization of std::locale
    // and pass for popular compilers and STL realizations
    BOOST_TEST(boost::lexical_cast<wide_char>(c_arr[0]) == wc_arr[0]);
    BOOST_TEST(boost::lexical_cast<wide_string>(c_arr) == wide_string(wc_arr));

    BOOST_TEST(boost::lexical_cast<wide_string>(sc_arr) == wide_string(wc_arr) );
    BOOST_TEST(boost::lexical_cast<wide_string>(uc_arr) == wide_string(wc_arr) );

    BOOST_TEST(boost::lexical_cast<wide_char>(uc_arr[0]) == wc_arr[0]);
    BOOST_TEST(boost::lexical_cast<wide_char>(sc_arr[0]) == wc_arr[0]);
}


void test_char_types_conversions_wchar_t()
{
#ifndef BOOST_LCAST_NO_WCHAR_T
    test_impl(L"Test array of chars");
    wchar_t c = boost::detail::lcast_char_constants<wchar_t>::zero;
    BOOST_TEST(L'0' == c);

    c = boost::detail::lcast_char_constants<wchar_t>::minus;
    BOOST_TEST(L'-' == c);

    c = boost::detail::lcast_char_constants<wchar_t>::plus;
    BOOST_TEST(L'+' == c);

    c = boost::detail::lcast_char_constants<wchar_t>::lowercase_e;
    BOOST_TEST(L'e' == c);

    c = boost::detail::lcast_char_constants<wchar_t>::capital_e;
    BOOST_TEST(L'E' == c);

    c = boost::detail::lcast_char_constants<wchar_t>::c_decimal_separator;
    BOOST_TEST(L'.' == c);
#endif

    BOOST_TEST(true);
}

void test_char_types_conversions_char16_t()
{
#if !defined(BOOST_NO_CXX11_CHAR16_T) && !defined(BOOST_NO_CXX11_UNICODE_LITERALS) && defined(BOOST_STL_SUPPORTS_NEW_UNICODE_LOCALES)
    test_impl(u"Test array of chars");
    char16_t c = boost::detail::lcast_char_constants<char16_t>::zero;
    BOOST_TEST(u'0' == c);

    c = boost::detail::lcast_char_constants<char16_t>::minus;
    BOOST_TEST(u'-' == c);

    c = boost::detail::lcast_char_constants<char16_t>::plus;
    BOOST_TEST(u'+' == c);

    c = boost::detail::lcast_char_constants<char16_t>::lowercase_e;
    BOOST_TEST(u'e' == c);

    c = boost::detail::lcast_char_constants<char16_t>::capital_e;
    BOOST_TEST(u'E' == c);

    c = boost::detail::lcast_char_constants<char16_t>::c_decimal_separator;
    BOOST_TEST(u'.' == c);
#endif

    BOOST_TEST(true);
}

void test_char_types_conversions_char32_t()
{
#if !defined(BOOST_NO_CXX11_CHAR32_T) && !defined(BOOST_NO_CXX11_UNICODE_LITERALS) && defined(BOOST_STL_SUPPORTS_NEW_UNICODE_LOCALES)
    test_impl(U"Test array of chars");
    char32_t c = boost::detail::lcast_char_constants<char32_t>::zero;
    BOOST_TEST(U'0' == c);

    c = boost::detail::lcast_char_constants<char32_t>::minus;
    BOOST_TEST(U'-' == c);

    c = boost::detail::lcast_char_constants<char32_t>::plus;
    BOOST_TEST(U'+' == c);

    c = boost::detail::lcast_char_constants<char32_t>::lowercase_e;
    BOOST_TEST(U'e' == c);

    c = boost::detail::lcast_char_constants<char32_t>::capital_e;
    BOOST_TEST(U'E', == c);

    c = boost::detail::lcast_char_constants<char32_t>::c_decimal_separator;
    BOOST_TEST(U'.' == c);
#endif

    BOOST_TEST(true);
}

int main()
{
    test_char_types_conversions_wchar_t();
    test_char_types_conversions_char16_t();
    test_char_types_conversions_char32_t();

    return boost::report_errors();
}
