//===----------------------------------------------------------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is dual licensed under the MIT and the University of Illinois Open
// Source Licenses. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//  Adaptation to Boost of the libcxx
//  Copyright 2010 Vicente J. Botet Escriba
//  Distributed under the Boost Software License, Version 1.0.
//  See http://www.boost.org/LICENSE_1_0.txt

#include <boost/ratio/ratio.hpp>
#include <cstdint>

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)

void test()
{
    {
    typedef boost::ratio<1, 1> R1;
    typedef boost::ratio<1, 1> R2;
    STATIC_ASSERT(!boost::ratio_not_equal<R1, R2>::value);
    }
    {
    typedef boost::ratio<INTMAX_MAX, 1> R1;
    typedef boost::ratio<INTMAX_MAX, 1> R2;
    STATIC_ASSERT(!boost::ratio_not_equal<R1, R2>::value);
    }
    {
    typedef boost::ratio<-INTMAX_MAX, 1> R1;
    typedef boost::ratio<-INTMAX_MAX, 1> R2;
    STATIC_ASSERT(!boost::ratio_not_equal<R1, R2>::value);
    }
    {
    typedef boost::ratio<1, INTMAX_MAX> R1;
    typedef boost::ratio<1, INTMAX_MAX> R2;
    STATIC_ASSERT(!boost::ratio_not_equal<R1, R2>::value);
    }
    {
    typedef boost::ratio<1, 1> R1;
    typedef boost::ratio<1, -1> R2;
    STATIC_ASSERT(boost::ratio_not_equal<R1, R2>::value);
    }
    {
    typedef boost::ratio<INTMAX_MAX, 1> R1;
    typedef boost::ratio<-INTMAX_MAX, 1> R2;
    STATIC_ASSERT(boost::ratio_not_equal<R1, R2>::value);
    }
    {
    typedef boost::ratio<-INTMAX_MAX, 1> R1;
    typedef boost::ratio<INTMAX_MAX, 1> R2;
    STATIC_ASSERT(boost::ratio_not_equal<R1, R2>::value);
    }
    {
    typedef boost::ratio<1, INTMAX_MAX> R1;
    typedef boost::ratio<1, -INTMAX_MAX> R2;
    STATIC_ASSERT(boost::ratio_not_equal<R1, R2>::value);
    }
}
