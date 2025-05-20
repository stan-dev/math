// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/charconv/limits.hpp>
#include <boost/charconv/detail/bit_layouts.hpp>

void test_odr_use( int const* );

template<typename T> void test()
{
    test_odr_use( &boost::charconv::limits<T>::max_chars10 );
    test_odr_use( &boost::charconv::limits<T>::max_chars );
}

void f2()
{
    test<char>();
    test<signed char>();
    test<unsigned char>();
    test<short>();
    test<unsigned short>();
    test<int>();
    test<unsigned>();
    test<long>();
    test<unsigned long>();
    test<long long>();
    test<unsigned long long>();

    test<float>();
    test<double>();
    #ifndef BOOST_CHARCONV_UNSUPPORTED_LONG_DOUBLE
    test<long double>();
    #endif

#ifdef BOOST_CHARCONV_HAS_INT128

    test<boost::int128_type>();
    test<boost::uint128_type>();

#endif
}
