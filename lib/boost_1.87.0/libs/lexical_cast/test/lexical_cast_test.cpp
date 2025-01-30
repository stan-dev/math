//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Terje Sletteb and Kevlin Henney, 2005.
//  Copyright Alexander Nasonov, 2006.
//  Copyright Antony Polukhin, 2011-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).
//
// Note: The unit test no longer compile on MSVC 6, but lexical_cast itself works for it.

//
// We need this #define before any #includes: otherwise msvc will emit warnings
// deep within std::string, resulting from our (perfectly legal) use of basic_string
// with a custom traits class:
//
#define _SCL_SECURE_NO_WARNINGS

#include <boost/lexical_cast.hpp>

#include <boost/cstdint.hpp>

#include <boost/core/lightweight_test.hpp>

#include <boost/type_traits/integral_promotion.hpp>
#include <string>
#include <vector>
#include <algorithm> // std::transform
#include <memory>

#if (defined(BOOST_HAS_LONG_LONG) || defined(BOOST_HAS_MS_INT64)) \
    && !(defined(BOOST_MSVC) && BOOST_MSVC < 1300)
#define LCAST_TEST_LONGLONG
#endif

#if defined(BOOST_NO_STRINGSTREAM) || defined(BOOST_NO_STD_WSTRING)
#define BOOST_LCAST_NO_WCHAR_T
#endif

#ifndef BOOST_TEST_CLOSE_FRACTION
// Naiive, but works for tests in this file 
#define BOOST_TEST_CLOSE_FRACTION(x, y, eps) BOOST_TEST(x - y + eps <= eps * 2)
#endif

template<class CharT>
struct my_traits : std::char_traits<CharT>
{
};

template<class CharT>
struct my_allocator : std::allocator<CharT>
{
    typedef std::allocator<CharT> base_t;

    my_allocator(){}
    template <class U> my_allocator(const my_allocator<U>& v) : base_t(v) {}

    template <class U> struct rebind { typedef my_allocator<U> other; };
};

using namespace boost;

void test_conversion_to_char()
{
    BOOST_TEST_EQ('A', lexical_cast<char>('A'));
    BOOST_TEST_EQ(' ', lexical_cast<char>(' '));
    BOOST_TEST_EQ('1', lexical_cast<char>(1));
    BOOST_TEST_EQ('0', lexical_cast<char>(0));
    BOOST_TEST_THROWS(lexical_cast<char>(123), bad_lexical_cast);
    BOOST_TEST_EQ('1', lexical_cast<char>(1.0));
    BOOST_TEST_EQ('1', lexical_cast<char>(true));
    BOOST_TEST_EQ('0', lexical_cast<char>(false));
    BOOST_TEST_EQ('A', lexical_cast<char>("A"));
    BOOST_TEST_EQ(' ', lexical_cast<char>(" "));
    BOOST_TEST_THROWS(lexical_cast<char>(""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<char>("Test"), bad_lexical_cast);
    BOOST_TEST_EQ('A', lexical_cast<char>(std::string("A")));
    BOOST_TEST_EQ(' ', lexical_cast<char>(std::string(" ")));
    BOOST_TEST_THROWS(
        lexical_cast<char>(std::string("")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<char>(std::string("Test")), bad_lexical_cast);
}

void test_conversion_to_int()
{
    BOOST_TEST_EQ(1, lexical_cast<int>('1'));
    BOOST_TEST_EQ(0, lexical_cast<int>('0'));
    BOOST_TEST_THROWS(lexical_cast<int>('A'), bad_lexical_cast);
    BOOST_TEST_EQ(1, lexical_cast<int>(1));
    BOOST_TEST_EQ(1, lexical_cast<int>(1.0));

    BOOST_TEST_EQ(
        (std::numeric_limits<int>::max)(),
        lexical_cast<int>((std::numeric_limits<int>::max)()));

    BOOST_TEST_EQ(
        (std::numeric_limits<int>::min)(),
        lexical_cast<int>((std::numeric_limits<int>::min)()));

    BOOST_TEST_THROWS(lexical_cast<int>(1.23), bad_lexical_cast);

    BOOST_TEST_THROWS(lexical_cast<int>(1e20), bad_lexical_cast);
    BOOST_TEST_EQ(1, lexical_cast<int>(true));
    BOOST_TEST_EQ(0, lexical_cast<int>(false));
    BOOST_TEST_EQ(123, lexical_cast<int>("123"));
    BOOST_TEST_THROWS(
        lexical_cast<int>(" 123"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<int>(""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<int>("Test"), bad_lexical_cast);
    BOOST_TEST_EQ(123, lexical_cast<int>("123"));
    BOOST_TEST_EQ(123, lexical_cast<int>(std::string("123")));
    BOOST_TEST_THROWS(
        lexical_cast<int>(std::string(" 123")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<int>(std::string("")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<int>(std::string("Test")), bad_lexical_cast);
}

void test_conversion_with_nonconst_char()
{
    std::vector<char> buffer;
    buffer.push_back('1');
    buffer.push_back('\0');
    BOOST_TEST_EQ(boost::lexical_cast<int>(&buffer[0]), 1);

    std::vector<unsigned char> buffer2;
    buffer2.push_back('1');
    buffer2.push_back('\0');
    BOOST_TEST_EQ(boost::lexical_cast<int>(&buffer2[0]), 1);

    std::vector<unsigned char> buffer3;
    buffer3.push_back('1');
    buffer3.push_back('\0');
    BOOST_TEST_EQ(boost::lexical_cast<int>(&buffer3[0]), 1);

#ifndef BOOST_LCAST_NO_WCHAR_T
    std::vector<wchar_t> buffer4;
    buffer4.push_back(L'1');
    buffer4.push_back(L'\0');
    BOOST_TEST_EQ(boost::lexical_cast<int>(&buffer4[0]), 1);
#endif
}

void test_conversion_to_double()
{
    BOOST_TEST_CLOSE_FRACTION(1.0, lexical_cast<double>('1'), (std::numeric_limits<double>::epsilon()));
    BOOST_TEST_THROWS(lexical_cast<double>('A'), bad_lexical_cast);
    BOOST_TEST_CLOSE_FRACTION(1.0, lexical_cast<double>(1), (std::numeric_limits<double>::epsilon()));
    BOOST_TEST_CLOSE_FRACTION(1.23, lexical_cast<double>(1.23), (std::numeric_limits<double>::epsilon()));
    BOOST_TEST_CLOSE_FRACTION(1.234567890, lexical_cast<double>(1.234567890), std::numeric_limits<double>::epsilon());
    BOOST_TEST_CLOSE_FRACTION(1.234567890, lexical_cast<double>("1.234567890"), std::numeric_limits<double>::epsilon());
    BOOST_TEST_CLOSE_FRACTION(1.0, lexical_cast<double>(true), (std::numeric_limits<double>::epsilon()));
    BOOST_TEST_CLOSE_FRACTION(0.0, lexical_cast<double>(false), (std::numeric_limits<double>::epsilon()));
    BOOST_TEST_CLOSE_FRACTION(1.23, lexical_cast<double>("1.23"), (std::numeric_limits<double>::epsilon()));
    BOOST_TEST_THROWS(lexical_cast<double>(""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<double>("Test"), bad_lexical_cast);
    BOOST_TEST_CLOSE_FRACTION(1.23, lexical_cast<double>(std::string("1.23")), (std::numeric_limits<double>::epsilon()));
    BOOST_TEST_THROWS(
        lexical_cast<double>(std::string("")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<double>(std::string("Test")), bad_lexical_cast);
}

void test_conversion_to_bool()
{
    BOOST_TEST_EQ(true, lexical_cast<bool>('1'));
    BOOST_TEST_EQ(false, lexical_cast<bool>('0'));
    BOOST_TEST_THROWS(lexical_cast<bool>('A'), bad_lexical_cast);
    BOOST_TEST_EQ(true, lexical_cast<bool>(1));
    BOOST_TEST_EQ(false, lexical_cast<bool>(0));
    BOOST_TEST_THROWS(lexical_cast<bool>(123), bad_lexical_cast);
    BOOST_TEST_EQ(true, lexical_cast<bool>(1.0));
    BOOST_TEST_THROWS(lexical_cast<bool>(-123), bad_lexical_cast);
    BOOST_TEST_EQ(false, lexical_cast<bool>(0.0));
    BOOST_TEST_THROWS(lexical_cast<bool>(1234), bad_lexical_cast);
#if !defined(_CRAYC)
    // Looks like a bug in CRAY compiler (throws bad_lexical_cast)
    // TODO: localize the bug and report it to developers.
    BOOST_TEST_EQ(true, lexical_cast<bool>(true));
    BOOST_TEST_EQ(false, lexical_cast<bool>(false));
#endif
    BOOST_TEST_EQ(true, lexical_cast<bool>("1"));
    BOOST_TEST_EQ(false, lexical_cast<bool>("0"));
    BOOST_TEST_THROWS(lexical_cast<bool>(""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("Test"), bad_lexical_cast);
    BOOST_TEST_EQ(true, lexical_cast<bool>("1"));
    BOOST_TEST_EQ(false, lexical_cast<bool>("0"));
    BOOST_TEST_EQ(true, lexical_cast<bool>(std::string("1")));
    BOOST_TEST_EQ(false, lexical_cast<bool>(std::string("0")));

    BOOST_TEST_THROWS(lexical_cast<bool>(1.0001L), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>(2), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>(2u), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>(-1), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>(-2), bad_lexical_cast);


    BOOST_TEST_THROWS(
        lexical_cast<bool>(std::string("")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<bool>(std::string("Test")), bad_lexical_cast);

    BOOST_TEST(lexical_cast<bool>("+1") == true);
    BOOST_TEST(lexical_cast<bool>("+0") == false);
    BOOST_TEST(lexical_cast<bool>("-0") == false);
    BOOST_TEST_THROWS(lexical_cast<bool>("--0"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("-+-0"), bad_lexical_cast);

    BOOST_TEST(lexical_cast<bool>("0") == false);
    BOOST_TEST(lexical_cast<bool>("1") == true);
    BOOST_TEST(lexical_cast<bool>("00") == false);
    BOOST_TEST(lexical_cast<bool>("00000000000") == false);
    BOOST_TEST(lexical_cast<bool>("000000000001") == true);
    BOOST_TEST(lexical_cast<bool>("+00") == false );
    BOOST_TEST(lexical_cast<bool>("-00") == false );
    BOOST_TEST(lexical_cast<bool>("+00000000001") == true );

    BOOST_TEST_THROWS(lexical_cast<bool>("020"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("00200"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("-00200"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("+00200"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("000000000002"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("-1"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("-0000000001"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("00000000011"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("001001"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("-00000000010"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>("-000000000100"), bad_lexical_cast);
}

void test_conversion_to_string()
{
    char buf[] = "hello";
    char* str = buf;
    BOOST_TEST_EQ(str, lexical_cast<std::string>(str));
    BOOST_TEST_EQ("A", lexical_cast<std::string>('A'));
    BOOST_TEST_EQ(" ", lexical_cast<std::string>(' '));
    BOOST_TEST_EQ("123", lexical_cast<std::string>(123));
    BOOST_TEST_EQ("1.23", lexical_cast<std::string>(1.23));
    BOOST_TEST_EQ("1.111111111", lexical_cast<std::string>(1.111111111));
    BOOST_TEST_EQ("1", lexical_cast<std::string>(true));
    BOOST_TEST_EQ("0", lexical_cast<std::string>(false));
    BOOST_TEST_EQ("Test", lexical_cast<std::string>("Test"));
    BOOST_TEST_EQ(" ", lexical_cast<std::string>(" "));
    BOOST_TEST_EQ("", lexical_cast<std::string>(""));
    BOOST_TEST_EQ("Test", lexical_cast<std::string>(std::string("Test")));
    BOOST_TEST_EQ(" ", lexical_cast<std::string>(std::string(" ")));
    BOOST_TEST_EQ("", lexical_cast<std::string>(std::string("")));
}

void test_conversion_from_to_wchar_t_alias()
{
    BOOST_TEST_EQ(123u, lexical_cast<unsigned short>("123"));
    BOOST_TEST_EQ(123u, lexical_cast<unsigned int>("123"));
    BOOST_TEST_EQ(123u, lexical_cast<unsigned long>("123"));
    BOOST_TEST_EQ(std::string("123"),
        lexical_cast<std::string>(static_cast<unsigned short>(123)));
    BOOST_TEST_EQ(std::string("123"), lexical_cast<std::string>(123u));
    BOOST_TEST_EQ(std::string("123"), lexical_cast<std::string>(123ul));
}

void test_conversion_from_wchar_t()
{
#ifndef BOOST_LCAST_NO_WCHAR_T
#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    BOOST_TEST_EQ(1, lexical_cast<int>(L'1'));
    BOOST_TEST_THROWS(lexical_cast<int>(L'A'), bad_lexical_cast);
#endif

    BOOST_TEST_EQ(123, lexical_cast<int>(L"123"));
    BOOST_TEST_THROWS(lexical_cast<int>(L""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<int>(L"Test"), bad_lexical_cast);

#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    BOOST_TEST_EQ(1.0, lexical_cast<double>(L'1'));
    BOOST_TEST_THROWS(lexical_cast<double>(L'A'), bad_lexical_cast);
#endif

    BOOST_TEST_EQ(1.23, lexical_cast<double>(L"1.23"));
    BOOST_TEST_THROWS(lexical_cast<double>(L""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<double>(L"Test"), bad_lexical_cast);

#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    BOOST_TEST_EQ(true, lexical_cast<bool>(L'1'));
    BOOST_TEST_EQ(false, lexical_cast<bool>(L'0'));
    BOOST_TEST_THROWS(lexical_cast<bool>(L'A'), bad_lexical_cast);
#endif
    BOOST_TEST_EQ(true, lexical_cast<bool>(L"1"));
    BOOST_TEST_EQ(false, lexical_cast<bool>(L"0"));
    BOOST_TEST_THROWS(lexical_cast<bool>(L""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<bool>(L"Test"), bad_lexical_cast);
#endif
}

void test_conversion_to_wchar_t()
{
#if !defined(BOOST_LCAST_NO_WCHAR_T) && !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    BOOST_TEST(L'1' == lexical_cast<wchar_t>(1));
    BOOST_TEST(L'0' == lexical_cast<wchar_t>(0));
    BOOST_TEST(L'1' == lexical_cast<wchar_t>('1'));
    BOOST_TEST(L'0' == lexical_cast<wchar_t>('0'));
    BOOST_TEST_THROWS(lexical_cast<wchar_t>(123), bad_lexical_cast);
    BOOST_TEST(L'1' == lexical_cast<wchar_t>(1.0));
    BOOST_TEST(L'0' == lexical_cast<wchar_t>(0.0));
    BOOST_TEST(L'1' == lexical_cast<wchar_t>(true));
    BOOST_TEST(L'0' == lexical_cast<wchar_t>(false));
    BOOST_TEST(L'A' == lexical_cast<wchar_t>(L'A'));
    BOOST_TEST(L' ' == lexical_cast<wchar_t>(L' '));
    BOOST_TEST(L'A' == lexical_cast<wchar_t>(L"A"));
    BOOST_TEST(L' ' == lexical_cast<wchar_t>(L" "));
    BOOST_TEST_THROWS(lexical_cast<wchar_t>(L""), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<wchar_t>(L"Test"), bad_lexical_cast);
    BOOST_TEST(L'A' == lexical_cast<wchar_t>(std::wstring(L"A")));
    BOOST_TEST(L' ' == lexical_cast<wchar_t>(std::wstring(L" ")));
    BOOST_TEST_THROWS(
        lexical_cast<wchar_t>(std::wstring(L"")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<wchar_t>(std::wstring(L"Test")), bad_lexical_cast);
#endif
    BOOST_TEST(true);
}

void test_conversion_from_wstring()
{
#ifndef BOOST_LCAST_NO_WCHAR_T
    BOOST_TEST_EQ(123, lexical_cast<int>(std::wstring(L"123")));
    BOOST_TEST_THROWS(
        lexical_cast<int>(std::wstring(L"")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<int>(std::wstring(L"Test")), bad_lexical_cast);

    BOOST_TEST_EQ(true, lexical_cast<bool>(std::wstring(L"1")));
    BOOST_TEST_EQ(false, lexical_cast<bool>(std::wstring(L"0")));
    BOOST_TEST_THROWS(
        lexical_cast<bool>(std::wstring(L"")), bad_lexical_cast);
    BOOST_TEST_THROWS(
        lexical_cast<bool>(std::wstring(L"Test")), bad_lexical_cast);
#endif
    BOOST_TEST(true);
}

void test_conversion_to_wstring()
{
#ifndef BOOST_LCAST_NO_WCHAR_T
    wchar_t buf[] = L"hello";
    wchar_t* str = buf;
    BOOST_TEST(str == lexical_cast<std::wstring>(str));
    BOOST_TEST(L"123" == lexical_cast<std::wstring>(123));
    BOOST_TEST(L"1.23" == lexical_cast<std::wstring>(1.23));
    BOOST_TEST(L"1" == lexical_cast<std::wstring>(true));
    BOOST_TEST(L"0" == lexical_cast<std::wstring>(false));
#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
    BOOST_TEST(L"A" == lexical_cast<std::wstring>(L'A'));
    BOOST_TEST(L" " == lexical_cast<std::wstring>(L' '));
    BOOST_TEST(L"A" == lexical_cast<std::wstring>('A'));
#endif
    BOOST_TEST(L"Test" == lexical_cast<std::wstring>(L"Test"));
    BOOST_TEST(L" " == lexical_cast<std::wstring>(L" "));
    BOOST_TEST(L"" == lexical_cast<std::wstring>(L""));
    BOOST_TEST(L"Test" == lexical_cast<std::wstring>(std::wstring(L"Test")));
    BOOST_TEST(L" " == lexical_cast<std::wstring>(std::wstring(L" ")));
    BOOST_TEST(L"" == lexical_cast<std::wstring>(std::wstring(L"")));
#endif
    BOOST_TEST(true);
}

void test_bad_lexical_cast()
{
    try
    {
        lexical_cast<int>(std::string("Test"));

        BOOST_TEST(false); // Exception expected
    }
    catch(const bad_lexical_cast &e)
    {
        BOOST_TEST(e.source_type() == typeid(std::string));
        BOOST_TEST(e.target_type() == typeid(int));
    }
}

void test_no_whitespace_stripping()
{
    BOOST_TEST_THROWS(lexical_cast<int>(" 123"), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<int>("123 "), bad_lexical_cast);
}

void test_traits()
{
    typedef std::basic_string<char, my_traits<char> > my_string;

    my_string const s("s");
    BOOST_TEST(boost::lexical_cast<char>(s) == s[0]);
    BOOST_TEST(boost::lexical_cast<my_string>(s) == s);
    BOOST_TEST(boost::lexical_cast<my_string>(-1) == "-1");
    BOOST_TEST(boost::lexical_cast<int>(my_string("42")) == 42);
    BOOST_TEST(boost::lexical_cast<double>(my_string("1.0")) == 1.0);
}

void test_wtraits()
{
    typedef std::basic_string<wchar_t, my_traits<wchar_t> > my_string;

    my_string const s(L"s");
    BOOST_TEST(boost::lexical_cast<wchar_t>(s) == s[0]);
    BOOST_TEST(boost::lexical_cast<my_string>(s) == s);
    BOOST_TEST(boost::lexical_cast<my_string>(-1) == L"-1");
    BOOST_TEST(boost::lexical_cast<int>(my_string(L"42")) == 42);
    BOOST_TEST(boost::lexical_cast<double>(my_string(L"1.0")) == 1.0);
}

void test_allocator()
{
// Following test cause compilation error on MSVC2012:
// (Reason: cannot convert from 'std::_Wrap_alloc<_Alloc>' to 'const my_allocator<CharT>')
//
// MSVC developer is notified about this issue
#if !defined(_MSC_VER) || (_MSC_VER < 1700)
    typedef std::basic_string< char
                             , std::char_traits<char>
                             , my_allocator<char>
                             > my_string;

    my_string s("s");
    BOOST_TEST(boost::lexical_cast<char>(s) == s[0]);
    BOOST_TEST(boost::lexical_cast<std::string>(s) == "s");
    BOOST_TEST(boost::lexical_cast<my_string>(s) == s);
    BOOST_TEST(boost::lexical_cast<my_string>(1) == "1");
    BOOST_TEST(boost::lexical_cast<my_string>("s") == s);
    BOOST_TEST(boost::lexical_cast<my_string>(std::string("s")) == s);
#endif
}

void test_wallocator()
{
// Following test cause compilation error on MSVC2012:
// (Reason: cannot convert from 'std::_Wrap_alloc<_Alloc>' to 'const my_allocator<CharT>')
//
// MSVC developer is notified about this issue
#if !defined(_MSC_VER) || (_MSC_VER < 1700)
    typedef std::basic_string< wchar_t
                             , std::char_traits<wchar_t>
                             , my_allocator<wchar_t>
                             > my_string;

    my_string s(L"s");
    BOOST_TEST(boost::lexical_cast<wchar_t>(s) == s[0]);
    BOOST_TEST(boost::lexical_cast<std::wstring>(s) == L"s");
    BOOST_TEST(boost::lexical_cast<my_string>(s) == s);
    BOOST_TEST(boost::lexical_cast<my_string>(1) == L"1");
    BOOST_TEST(boost::lexical_cast<my_string>(L"s") == s);
    BOOST_TEST(boost::lexical_cast<my_string>(std::wstring(L"s")) == s);
#endif
}


void test_char_types_conversions()
{
    const char c_arr[]            = "Test array of chars";
    const unsigned char uc_arr[]  = "Test array of chars";
    const signed char sc_arr[]    = "Test array of chars";

    BOOST_TEST(boost::lexical_cast<std::string>(c_arr) == std::string(c_arr));
    BOOST_TEST(boost::lexical_cast<std::string>(uc_arr) == std::string(c_arr));
    BOOST_TEST(boost::lexical_cast<std::string>(sc_arr) == std::string(c_arr));

    BOOST_TEST(boost::lexical_cast<char>(c_arr[0]) == c_arr[0]);
    BOOST_TEST(boost::lexical_cast<char>(uc_arr[0]) == c_arr[0]);
    BOOST_TEST(boost::lexical_cast<char>(sc_arr[0]) == c_arr[0]);

    BOOST_TEST(boost::lexical_cast<unsigned char>(c_arr[0]) == uc_arr[0]);
    BOOST_TEST(boost::lexical_cast<unsigned char>(uc_arr[0]) == uc_arr[0]);
    BOOST_TEST(boost::lexical_cast<unsigned char>(sc_arr[0]) == uc_arr[0]);

    BOOST_TEST(boost::lexical_cast<signed char>(c_arr[0]) == sc_arr[0]);
    BOOST_TEST(boost::lexical_cast<signed char>(uc_arr[0]) == sc_arr[0]);
    BOOST_TEST(boost::lexical_cast<signed char>(sc_arr[0]) == sc_arr[0]);

#ifndef BOOST_LCAST_NO_WCHAR_T
    const wchar_t wc_arr[]=L"Test array of chars";

    BOOST_TEST(boost::lexical_cast<std::wstring>(wc_arr) == std::wstring(wc_arr));
    BOOST_TEST(boost::lexical_cast<wchar_t>(wc_arr[0]) == wc_arr[0]);

#endif
}



struct foo_operators_test
{
  foo_operators_test() : f(2) {}
  int f;
};

template <typename OStream>
OStream& operator<<(OStream& ostr, const foo_operators_test& foo)
{
  ostr << foo.f;
  return ostr;
}

template <typename IStream>
IStream& operator>>(IStream& istr, foo_operators_test& foo)
{
  istr >> foo.f;
  return istr;
}

void operators_overload_test()
{
    foo_operators_test foo;
    BOOST_TEST_EQ(boost::lexical_cast<std::string>(foo), "2");
    BOOST_TEST_EQ((boost::lexical_cast<foo_operators_test>("2")).f, 2);

    // Must compile
    (void)boost::lexical_cast<foo_operators_test>(foo);
}


void test_char16_conversions()
{
// There's no std::ctype<char16_t> in Xcode_15.0.1
#if !defined(BOOST_NO_CXX11_CHAR16_T) && !defined(BOOST_NO_CXX11_UNICODE_LITERALS) && !defined(__APPLE__)
    BOOST_TEST(u"100" == lexical_cast<std::u16string>(u"100"));
    BOOST_TEST(u"1" == lexical_cast<std::u16string>(u'1'));
#endif
}

void test_char32_conversions()
{
// There's no std::ctype<char32_t> in Xcode_15.0.1
#if !defined(BOOST_NO_CXX11_CHAR32_T) && !defined(BOOST_NO_CXX11_UNICODE_LITERALS) && !defined(__APPLE__)
    BOOST_TEST(U"100" == lexical_cast<std::u32string>(U"100"));
    BOOST_TEST(U"1" == lexical_cast<std::u32string>(U'1'));
#endif
}

void test_getting_pointer_to_function()
{
    // Just checking that &lexical_cast<To, From> is not ambiguous
    typedef char char_arr[4];
    typedef int(*f1)(const char_arr&);
    f1 p1 = &boost::lexical_cast<int, char_arr>;
    BOOST_TEST(p1);

    typedef int(*f2)(const std::string&);
    f2 p2 = &boost::lexical_cast<int, std::string>;
    BOOST_TEST(p2);

    typedef std::string(*f3)(const int&);
    f3 p3 = &boost::lexical_cast<std::string, int>;
    BOOST_TEST(p3);

    std::vector<int> values;
    std::vector<std::string> ret;
    std::transform(values.begin(), values.end(), ret.begin(), boost::lexical_cast<std::string, int>);
}

int main()
{
    test_conversion_to_char();
    test_conversion_to_int();
    test_conversion_to_double();
    test_conversion_to_bool();
    test_conversion_from_to_wchar_t_alias();
    test_conversion_to_string();
    test_conversion_with_nonconst_char();
#ifndef BOOST_LCAST_NO_WCHAR_T
    test_conversion_from_wchar_t();
    test_conversion_to_wchar_t();
    test_conversion_from_wstring();
    test_conversion_to_wstring();
#endif
    test_bad_lexical_cast();
    test_no_whitespace_stripping();

    test_traits();
    test_wtraits();
    test_allocator();
    test_wallocator();

    test_char_types_conversions();
    operators_overload_test();
    test_char16_conversions();
    test_char32_conversions();
    test_getting_pointer_to_function();

    return boost::report_errors();
}


