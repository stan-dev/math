// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/ratio/ratio.hpp>

int main()
{
    typedef boost::ratio<5, 2> R1;
    typedef boost::ratio<2, 7> R2;
    typedef boost::ratio_multiply<R1, R2>::type R3;

    return R3::num == 5 && R3::den == 7? 0: 1;
}
