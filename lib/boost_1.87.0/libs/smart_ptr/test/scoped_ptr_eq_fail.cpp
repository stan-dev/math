//
// scoped_ptr_eq_fail.cpp - a negative test for "p == q"
//
// Copyright (c) 2008 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/scoped_ptr.hpp>

int main()
{
    boost::scoped_ptr<int> p, q;
    p == q; // must fail
    return 0;
}
