//
//  shared_ptr_pv_fail.cpp - a negative test for converting a shared_ptr to void*
//
//  Copyright 2007 Peter Dimov
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt
//

#include <boost/shared_ptr.hpp>

void f( void* )
{
}

int main()
{
    boost::shared_ptr<int> p;
    f( p ); // must fail
    return 0;
}
