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

template <typename R>
struct numerator;

template <std::intmax_t N, std::intmax_t D>
struct numerator<boost::ratio<N,D> > {
    static const std::intmax_t value = N;
};

STATIC_ASSERT(numerator<boost::ratio_add<boost::ratio<1,2>,boost::ratio<1,3> > >::value == 1);
