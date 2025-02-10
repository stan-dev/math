// Copyright 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 8
# pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <vector>
#include <functional> // to catch msvc-14.1 conflicts

template<class T> void test()
{
    typedef std::vector<T> list;
    typedef boost::hash<list> hash;

    int const N = 32;

    std::size_t h[ N ];

    list v;

    for( int i = 0; i < N; ++i )
    {
        h[ i ] = hash()( v );

        BOOST_TEST_EQ( h[ i ], hash()( v ) );

        for( int j = 0; j < i; ++j )
        {
            BOOST_TEST_NE( h[ j ], h[ i ] );
        }

        v.push_back( T() );
    }
}

int main()
{
    test<int>();
    test<float>();
    test<double>();
    test< std::vector<int> >();

    return boost::report_errors();
}
