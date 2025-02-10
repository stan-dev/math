//
//  auto_ptr_lv_fail.cpp - a negative test for converting an auto_ptr to shared_ptr
//
//  Copyright 2009 Peter Dimov
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt
//

#include <boost/shared_ptr.hpp>
#include <memory>

void f( boost::shared_ptr<int> )
{
}

int main()
{
    std::auto_ptr<int> p;
    f( p ); // must fail
    return 0;
}
