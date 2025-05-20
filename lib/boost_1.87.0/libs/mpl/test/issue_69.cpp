// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0. 
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mpl/integral_c.hpp>

enum test { zero, one };

int main()
{
    boost::mpl::integral_c<test, zero> m;
    return m.value;
}
