// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <vector>
#include <deque>
#include <list>
#if !defined(BOOST_NO_CXX11_HDR_FORWARD_LIST)
# include <forward_list>
#endif

template<class T> std::size_t hv( T const& t )
{
    return boost::hash<T>()( t );
}

template<class T> void test()
{
    for( std::size_t i = 0; i < 8; ++i )
    {
        std::vector<T> v( i );
        std::size_t h0 = hv( v );

        std::deque<T> d( v.begin(), v.end() );
        BOOST_TEST_EQ( h0, hv( d ) );

        std::list<T> l( v.begin(), v.end() );
        BOOST_TEST_EQ( h0, hv( l ) );

#if !defined(BOOST_NO_CXX11_HDR_FORWARD_LIST)

        std::forward_list<T> f( v.begin(), v.end() );
        BOOST_TEST_EQ( h0, hv( f ) );

#endif
    }
}

int main()
{
    test<char>();
    test<unsigned char>();
    test<signed char>();
    test<int>();
    test<float>();
    test<double>();

#if defined(__cpp_char8_t) && __cpp_char8_t >= 201811L
    test<char8_t>();
#endif

#if defined(__cpp_lib_byte) && __cpp_lib_byte >= 201603L
    test<std::byte>();
#endif

    return boost::report_errors();
}
