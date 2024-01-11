// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/ratio/ratio.hpp>

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)

void test()
{
    {
    typedef boost::ratio<1, 1> R1;
    typedef boost::ratio<1, 1> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 1);
    }
    {
    typedef boost::ratio<1, 2> R1;
    typedef boost::ratio<1, 1> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 2);
    }
    {
    typedef boost::ratio<-1, 2> R1;
    typedef boost::ratio<1, 1> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 2);
    }
    {
    typedef boost::ratio<1, -2> R1;
    typedef boost::ratio<1, 1> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 2);
    }
    {
    typedef boost::ratio<1, 2> R1;
    typedef boost::ratio<-1, 1> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 2);
    }
    {
    typedef boost::ratio<1, 2> R1;
    typedef boost::ratio<1, -1> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 2);
    }
    {
    typedef boost::ratio<1, 2> R1;
    typedef boost::ratio<1, 3> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 6);
    }
    {
    typedef boost::ratio<7, 13> R1;
    typedef boost::ratio<11, 13> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 1 && R::den == 13);
    }
    {
    typedef boost::ratio<21, 55> R1;
    typedef boost::ratio<14, 25> R2;
    typedef boost::ratio_gcd<R1, R2> R;
    STATIC_ASSERT(R::num == 7 && R::den == 275);
    }
}
