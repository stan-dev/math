// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>
#include <cstddef>

#if defined(BOOST_NO_CXX11_NULLPTR)

BOOST_PRAGMA_MESSAGE( "Test skipped, BOOST_NO_CXX11_NULLPTR is defined" )
int main() {}

#else

template<class T> std::size_t hv( T const& x )
{
    return boost::hash<T>()( x );
}

int main()
{
    {
        BOOST_TEST_EQ( hv((void*)0), hv(nullptr) );
    }

    {
        int x = 0;

        BOOST_TEST_EQ( hv((int*)0), hv(nullptr) );
        BOOST_TEST_NE( hv(&x), hv(nullptr) );

        (void)x;
    }

    return boost::report_errors();
}

#endif
