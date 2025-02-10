#include <boost/config.hpp>

#if defined(BOOST_MSVC)
#pragma warning(disable: 4786)  // identifier truncated in debug info
#pragma warning(disable: 4710)  // function not inlined
#pragma warning(disable: 4711)  // function selected for automatic inline expansion
#pragma warning(disable: 4514)  // unreferenced inline removed
#endif

//
// mem_fn_dm_test.cpp - data members
//
// Copyright 2005, 2024 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mem_fn.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/core/ref.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config/workaround.hpp>

struct X
{
    int m;
};

int main()
{
    X x = { 0 };

    boost::mem_fn( &X::m )( x ) = 401;

    BOOST_TEST( x.m == 401 );
    BOOST_TEST( boost::mem_fn( &X::m )( x ) == 401 );

    boost::mem_fn( &X::m )( &x ) = 502;

    BOOST_TEST( x.m == 502 );
    BOOST_TEST( boost::mem_fn( &X::m )( &x ) == 502 );

    X * px = &x;

    boost::mem_fn( &X::m )( px ) = 603;

    BOOST_TEST( x.m == 603 );
    BOOST_TEST( boost::mem_fn( &X::m )( px ) == 603 );

    X const & cx = x;
    X const * pcx = &x;

    BOOST_TEST( boost::mem_fn( &X::m )( cx ) == 603 );
    BOOST_TEST( boost::mem_fn( &X::m )( pcx ) == 603 );

    BOOST_TEST_EQ( boost::mem_fn( &X::m )( boost::ref( x ) ), 603 );
    BOOST_TEST_EQ( boost::mem_fn( &X::m )( boost::cref( x ) ), 603 );

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1900)

    boost::mem_fn( &X::m )( boost::ref( x ) ) = 704;
    BOOST_TEST_EQ( x.m, 704 );

#endif

    boost::shared_ptr<X> sp( new X() );
    boost::shared_ptr<X const> csp( sp );

    sp->m = 805;

    BOOST_TEST_EQ( boost::mem_fn( &X::m )( sp ), 805 );
    BOOST_TEST_EQ( boost::mem_fn( &X::m )( csp ), 805 );

#if !BOOST_WORKAROUND(BOOST_MSVC, < 1900)

    boost::mem_fn( &X::m )( sp ) = 906;
    BOOST_TEST_EQ( sp->m, 906 );

#endif

    return boost::report_errors();
}
