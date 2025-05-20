// Copyright 2017, 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/typeof/typeof.hpp>

int main()
{
    BOOST_AUTO(x, 5);
    return x == 5? 0: 1;
}
