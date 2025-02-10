//
//  shared_ptr_assign_fail.cpp - a negative test for shared_ptr assignment
//
//  Copyright (c) 2002 Peter Dimov and Multi Media Ltd.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/shared_ptr.hpp>

int main()
{
    boost::shared_ptr<int> p;
    p = new int(42); // assignment must fail
    return 0;
}
