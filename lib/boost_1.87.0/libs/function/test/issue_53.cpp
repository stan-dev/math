// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// https://github.com/boostorg/function/issues/53

#include <boost/bind/apply.hpp>
#include <boost/bind/bind.hpp>
#include <boost/function.hpp>

int TestArg( int, double )
{
    return 0;
}


void f()
{
    boost::function<int(int)> fn = boost::bind( &TestArg, boost::placeholders::_1, 1.0 );
}
