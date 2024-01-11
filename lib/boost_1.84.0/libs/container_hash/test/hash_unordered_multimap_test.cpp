// Copyright 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 8
# pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>

#if defined(BOOST_NO_CXX11_HDR_UNORDERED_MAP)

BOOST_PRAGMA_MESSAGE( "Test skipped, BOOST_NO_CXX11_HDR_UNORDERED_MAP is defined" )
int main() {}

#elif defined(BOOST_NO_CXX11_VARIADIC_TEMPLATES) && defined(_CPPLIB_VER) && _CPPLIB_VER >= 520

BOOST_PRAGMA_MESSAGE( "Test skipped, _CPPLIB_VER >= 520 and BOOST_NO_CXX11_VARIADIC_TEMPLATES is defined" )
int main() {}

#else

#include <unordered_map>

template<class T> void test()
{
    typedef std::unordered_multimap<T, T> map;
    typedef boost::hash<map> hash;

    int const N = 32;

    std::size_t h[ N ];

    map v;

    for( int i = 0; i < N; ++i )
    {
        h[ i ] = hash()( v );

        BOOST_TEST_EQ( h[ i ], hash()( v ) );

        for( int j = 0; j < i; ++j )
        {
            BOOST_TEST_NE( h[ j ], h[ i ] );
        }

        v.insert( std::pair<T const, T>() );
    }
}

int main()
{
    test<int>();
    test<float>();
    test<double>();

    return boost::report_errors();
}

#endif
