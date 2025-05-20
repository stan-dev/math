//
//  shared_ptr_delete_fail.cpp - a negative test for "delete sp;"
//
//  Copyright 2005 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/shared_ptr.hpp>

int main()
{
    boost::shared_ptr<int> p;
    delete p; // must fail
    return 0;
}
