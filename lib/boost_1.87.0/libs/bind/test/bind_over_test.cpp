//
// bind_over_test.cpp - overloaded functions
//
// Copyright 2017, 2024 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/bind/bind.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::placeholders;

//

int f()
{
    return 17041;
}

int f( int x1 )
{
    return x1;
}

int f( int x1, int x2 )
{
    return x1+x2;
}

int f( int x1, int x2, int x3 )
{
    return x1+x2+x3;
}

int f( int x1, int x2, int x3, int x4 )
{
    return x1+x2+x3+x4;
}

int f( int x1, int x2, int x3, int x4, int x5 )
{
    return x1+x2+x3+x4+x5;
}

int f( int x1, int x2, int x3, int x4, int x5, int x6 )
{
    return x1+x2+x3+x4+x5+x6;
}

int f( int x1, int x2, int x3, int x4, int x5, int x6, int x7 )
{
    return x1+x2+x3+x4+x5+x6+x7;
}

int f( int x1, int x2, int x3, int x4, int x5, int x6, int x7, int x8 )
{
    return x1+x2+x3+x4+x5+x6+x7+x8;
}

int f( int x1, int x2, int x3, int x4, int x5, int x6, int x7, int x8, int x9 )
{
    return x1+x2+x3+x4+x5+x6+x7+x8+x9;
}

void overloaded_function_test()
{
    BOOST_TEST_EQ( boost::bind( f )(), 17041 );
    BOOST_TEST_EQ( boost::bind( f, 1 )(), 1 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2 )(), 1+2 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2, 3 )(), 1+2+3 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2, 3, 4 )(), 1+2+3+4 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2, 3, 4, 5 )(), 1+2+3+4+5 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2, 3, 4, 5, 6 )(), 1+2+3+4+5+6 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2, 3, 4, 5, 6, 7 )(), 1+2+3+4+5+6+7 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2, 3, 4, 5, 6, 7, 8 )(), 1+2+3+4+5+6+7+8 );
    BOOST_TEST_EQ( boost::bind( f, 1, 2, 3, 4, 5, 6, 7, 8, 9 )(), 1+2+3+4+5+6+7+8+9 );
}

int main()
{
    overloaded_function_test();
    return boost::report_errors();
}
