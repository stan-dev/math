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

#include <boost/core/cmath.hpp>
#include <boost/type_traits/is_same.hpp>

#include <boost/core/lightweight_test.hpp>

#if defined(BOOST_NO_STRINGSTREAM) || defined(BOOST_NO_STD_WSTRING)
#define BOOST_LCAST_NO_WCHAR_T
#endif

using namespace boost;

template <class T>
bool is_pos_inf(T value)
{
    return (boost::core::isinf)(value) && !(boost::core::signbit)(value);
}

template <class T>
bool is_neg_inf(T value)
{
    return (boost::core::isinf)(value) && (boost::core::signbit)(value);
}

template <class T>
bool is_pos_nan(T value)
{
    return (boost::core::isnan)(value) && !(boost::core::signbit)(value);
}

template <class T>
bool is_neg_nan(T value)
{
    /* There is some strange behaviour on Itanium platform with -nan nuber for long double.
    * It is a IA64 feature, or it is a boost::math feature, not a lexical_cast bug */
#if defined(__ia64__) || defined(_M_IA64)
    return (boost::core::isnan)(value)
            && ( boost::is_same<T, long double >::value || (boost::core::signbit)(value) );
#else
    return (boost::core::isnan)(value) && (boost::core::signbit)(value);
#endif
}

template <class T>
void test_inf_nan_templated()
{
    typedef T test_t;

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("inf") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("INF") ) );

    BOOST_TEST( is_neg_inf( lexical_cast<test_t>("-inf") ) );
    BOOST_TEST( is_neg_inf( lexical_cast<test_t>("-INF") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("+inf") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("+INF") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("infinity") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("INFINITY") ) );

    BOOST_TEST( is_neg_inf( lexical_cast<test_t>("-infinity") ) );
    BOOST_TEST( is_neg_inf( lexical_cast<test_t>("-INFINITY") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("+infinity") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("+INFINITY") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("iNfiNity") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>("INfinity") ) );

    BOOST_TEST( is_neg_inf( lexical_cast<test_t>("-inFINITY") ) );
    BOOST_TEST( is_neg_inf( lexical_cast<test_t>("-INFINITY") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("nan") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("NAN") ) );

    BOOST_TEST( is_neg_nan( lexical_cast<test_t>("-nan") ) );
    BOOST_TEST( is_neg_nan( lexical_cast<test_t>("-NAN") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("+nan") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("+NAN") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("nAn") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("NaN") ) );

    BOOST_TEST( is_neg_nan( lexical_cast<test_t>("-nAn") ) );
    BOOST_TEST( is_neg_nan( lexical_cast<test_t>("-NaN") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("+Nan") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("+nAN") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("nan()") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>("NAN(some string)") ) );
    BOOST_TEST_THROWS( lexical_cast<test_t>("NAN(some string"), bad_lexical_cast );

    BOOST_TEST(lexical_cast<std::string>( (boost::core::copysign)(std::numeric_limits<test_t >::infinity(), static_cast<test_t>(-1.0)))
                == "-inf" );
    BOOST_TEST(lexical_cast<std::string>( std::numeric_limits<test_t >::infinity()) == "inf" );
    BOOST_TEST(lexical_cast<std::string>( std::numeric_limits<test_t >::quiet_NaN()) == "nan" );
#if !defined(__ia64__) && !defined(_M_IA64)
    BOOST_TEST(lexical_cast<std::string>(
                (boost::core::copysign)(std::numeric_limits<test_t >::quiet_NaN(), static_cast<test_t>(-1.0)))
                == "-nan" );
#endif

#ifndef  BOOST_LCAST_NO_WCHAR_T
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"inf") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"INF") ) );

    BOOST_TEST( is_neg_inf( lexical_cast<test_t>(L"-inf") ) );
    BOOST_TEST( is_neg_inf( lexical_cast<test_t>(L"-INF") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"+inf") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"+INF") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"infinity") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"INFINITY") ) );

    BOOST_TEST( is_neg_inf( lexical_cast<test_t>(L"-infinity") ) );
    BOOST_TEST( is_neg_inf( lexical_cast<test_t>(L"-INFINITY") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"+infinity") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"+INFINITY") ) );

    BOOST_TEST( is_neg_inf( lexical_cast<test_t>(L"-infINIty") ) );
    BOOST_TEST( is_neg_inf( lexical_cast<test_t>(L"-INFiniTY") ) );

    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"+inFINIty") ) );
    BOOST_TEST( is_pos_inf( lexical_cast<test_t>(L"+INfinITY") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>(L"nan") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>(L"NAN") ) );

    BOOST_TEST( is_neg_nan( lexical_cast<test_t>(L"-nan") ) );
    BOOST_TEST( is_neg_nan( lexical_cast<test_t>(L"-NAN") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>(L"+nan") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>(L"+NAN") ) );

    BOOST_TEST( is_pos_nan( lexical_cast<test_t>(L"nan()") ) );
    BOOST_TEST( is_pos_nan( lexical_cast<test_t>(L"NAN(some string)") ) );
    BOOST_TEST_THROWS( lexical_cast<test_t>(L"NAN(some string"), bad_lexical_cast );

    BOOST_TEST(lexical_cast<std::wstring>( (boost::core::copysign)(std::numeric_limits<test_t >::infinity(), static_cast<test_t>(-1.0)))
                == L"-inf" );
    BOOST_TEST(lexical_cast<std::wstring>( std::numeric_limits<test_t >::infinity()) == L"inf" );
    BOOST_TEST(lexical_cast<std::wstring>( std::numeric_limits<test_t >::quiet_NaN()) == L"nan" );
#if !defined(__ia64__) && !defined(_M_IA64)
    BOOST_TEST(lexical_cast<std::wstring>(
                (boost::core::copysign)(std::numeric_limits<test_t >::quiet_NaN(), static_cast<test_t>(-1.0)))
                == L"-nan" );
#endif

#endif
}

void test_inf_nan_float()
{
    test_inf_nan_templated<float >();

    BOOST_TEST(is_pos_nan(lexical_cast<float>(std::numeric_limits<float>::quiet_NaN())));
    BOOST_TEST(is_neg_nan(lexical_cast<float>(-std::numeric_limits<float>::quiet_NaN())));
    BOOST_TEST(is_pos_inf(lexical_cast<float>(std::numeric_limits<float>::infinity())));
    BOOST_TEST(is_neg_inf(lexical_cast<float>(-std::numeric_limits<float>::infinity())));

    BOOST_TEST(is_pos_nan(lexical_cast<float>(std::numeric_limits<double>::quiet_NaN())));
    BOOST_TEST(is_neg_nan(lexical_cast<float>(-std::numeric_limits<double>::quiet_NaN())));
    BOOST_TEST(is_pos_inf(lexical_cast<float>(std::numeric_limits<double>::infinity())));
    BOOST_TEST(is_neg_inf(lexical_cast<float>(-std::numeric_limits<double>::infinity())));
}

void test_inf_nan_double()
{
    test_inf_nan_templated<double >();

    BOOST_TEST(is_pos_nan(lexical_cast<double>(std::numeric_limits<float>::quiet_NaN())));
    BOOST_TEST(is_neg_nan(lexical_cast<double>(-std::numeric_limits<float>::quiet_NaN())));
    BOOST_TEST(is_pos_inf(lexical_cast<double>(std::numeric_limits<float>::infinity())));
    BOOST_TEST(is_neg_inf(lexical_cast<double>(-std::numeric_limits<float>::infinity())));

    BOOST_TEST(is_pos_nan(lexical_cast<double>(std::numeric_limits<double>::quiet_NaN())));
    BOOST_TEST(is_neg_nan(lexical_cast<double>(-std::numeric_limits<double>::quiet_NaN())));
    BOOST_TEST(is_pos_inf(lexical_cast<double>(std::numeric_limits<double>::infinity())));
    BOOST_TEST(is_neg_inf(lexical_cast<double>(-std::numeric_limits<double>::infinity())));

    BOOST_TEST(is_pos_nan(lexical_cast<double>(std::numeric_limits<long double>::quiet_NaN())));
    BOOST_TEST(is_neg_nan(lexical_cast<double>(-std::numeric_limits<long double>::quiet_NaN())));
    BOOST_TEST(is_pos_inf(lexical_cast<double>(std::numeric_limits<long double>::infinity())));
    BOOST_TEST(is_neg_inf(lexical_cast<double>(-std::numeric_limits<long double>::infinity())));
}

void test_inf_nan_long_double()
{
// We do not run tests on compilers with bugs
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_inf_nan_templated<long double >();
#endif
    BOOST_TEST(true);
}

int main()
{
    test_inf_nan_float();
    test_inf_nan_double();
    test_inf_nan_long_double();

    return boost::report_errors();
}
