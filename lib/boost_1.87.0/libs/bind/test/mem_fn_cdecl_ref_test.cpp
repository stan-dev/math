#include <boost/config.hpp>

#ifndef _MSC_VER

int main()
{
}

#else

//
//  mem_fn_cdecl_ref_test.cpp - reference_wrapper, cdecl
//
//  Copyright 2009, 2024 Peter Dimov
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt
//

#define BOOST_MEM_FN_ENABLE_CDECL

#include <boost/mem_fn.hpp>
#include <boost/ref.hpp>
#include <boost/core/lightweight_test.hpp>

struct X
{
    int __cdecl f0() { return 1; }
    int __cdecl g0() const { return -1; }

    int __cdecl f1( int x1 ) { return x1; }
    int __cdecl g1( int x1 ) const { return -x1; }

    int __cdecl f2( int x1, int x2 ) { return x1+x2; }
    int __cdecl g2( int x1, int x2 ) const { return -x1-x2; }

    int __cdecl f3( int x1, int x2, int x3 ) { return x1+x2+x3; }
    int __cdecl g3( int x1, int x2, int x3 ) const { return -x1-x2-x3; }
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

#endif
