//
//  mem_fn_ref_test.cpp - reference_wrapper
//
//  Copyright 2009, 2024 Peter Dimov
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt
//

#if defined(__GNUC__) && defined(__SANITIZE_ADDRESS__)
# pragma GCC diagnostic ignored "-Warray-bounds"
#endif

#include <boost/mem_fn.hpp>
#include <boost/ref.hpp>
#include <boost/core/lightweight_test.hpp>

struct X
{
    int f0() { return 1; }
    int g0() const { return -1; }

    int f1( int x1 ) { return x1; }
    int g1( int x1 ) const { return -x1; }

    int f2( int x1, int x2 ) { return x1+x2; }
    int g2( int x1, int x2 ) const { return -x1-x2; }

    int f3( int x1, int x2, int x3 ) { return x1+x2+x3; }
    int g3( int x1, int x2, int x3 ) const { return -x1-x2-x3; }
};

int main()
{
    X x;

    BOOST_TEST_EQ( boost::mem_fn( &X::f0 )( boost::ref( x ) ), 1 );
    BOOST_TEST_EQ( boost::mem_fn( &X::g0 )( boost::cref( x ) ), -1 );

    BOOST_TEST_EQ( boost::mem_fn( &X::f1 )( boost::ref( x ), 1 ), 1 );
    BOOST_TEST_EQ( boost::mem_fn( &X::g1 )( boost::cref( x ), 1 ), -1 );

    BOOST_TEST_EQ( boost::mem_fn( &X::f2 )( boost::ref( x ), 1, 2 ), 1+2 );
    BOOST_TEST_EQ( boost::mem_fn( &X::g2 )( boost::cref( x ), 1, 2 ), -1-2 );

    BOOST_TEST_EQ( boost::mem_fn( &X::f3 )( boost::ref( x ), 1, 2, 3 ), 1+2+3 );
    BOOST_TEST_EQ( boost::mem_fn( &X::g3 )( boost::cref( x ), 1, 2, 3 ), -1-2-3 );

    return boost::report_errors();
}
