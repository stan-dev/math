// Copyright 2005-2009 Daniel James.
// Copyright 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/type_name.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/config.hpp>

#if defined(BOOST_MSVC)
#pragma warning(disable: 4127) // conditional expression is constant
#endif

// This test checks that small numbers hash to themselves even if
// their type is wider than size_t

template<class T>
void identity_test()
{
    if( boost::is_signed<T>::value )
    {
        for( int i = -128; i <= 127; ++i )
        {
            BOOST_TEST_EQ( boost::hash<T>()( static_cast<T>( i ) ), static_cast<std::size_t>( i ) );
        }
    }
    else
    {
        for( int i = 0; i <= 255; ++i )
        {
            BOOST_TEST_EQ( boost::hash<T>()( static_cast<T>( i ) ), static_cast<std::size_t>( i ) );
        }
    }
}

#define TEST(type) std::cerr << "Testing: " #type " (" << boost::core::type_name<type>() << ")\n"; identity_test<type>();

int main()
{
    TEST(char)
    TEST(signed char)
    TEST(unsigned char)
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
    TEST(wchar_t)
#endif
#ifndef BOOST_NO_CXX11_CHAR16_T
    TEST(char16_t)
#endif
#ifndef BOOST_NO_CXX11_CHAR32_T
    TEST(char32_t)
#endif
    TEST(short)
    TEST(unsigned short)
    TEST(int)
    TEST(unsigned int)
    TEST(long)
    TEST(unsigned long)

#if !defined(BOOST_NO_LONG_LONG)
    TEST(boost::long_long_type)
    TEST(boost::ulong_long_type)
#endif

#if defined(BOOST_HAS_INT128)
    TEST(boost::int128_type)
    TEST(boost::uint128_type)
#endif

    TEST(std::size_t)
    TEST(std::ptrdiff_t)

    return boost::report_errors();
}
