//  Copyright 2010 Vicente J. Botet Escriba
//  Distributed under the Boost Software License, Version 1.0.
//  See http://www.boost.org/LICENSE_1_0.txt

// test ratio typedef's

#include <boost/ratio/ratio.hpp>

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)

STATIC_ASSERT(boost::atto::num == 1 && boost::atto::den == 1000000000000000000ULL);
STATIC_ASSERT(boost::femto::num == 1 && boost::femto::den == 1000000000000000ULL);
STATIC_ASSERT(boost::pico::num == 1 && boost::pico::den == 1000000000000ULL);
STATIC_ASSERT(boost::nano::num == 1 && boost::nano::den == 1000000000ULL);
STATIC_ASSERT(boost::micro::num == 1 && boost::micro::den == 1000000ULL);
STATIC_ASSERT(boost::milli::num == 1 && boost::milli::den == 1000ULL);
STATIC_ASSERT(boost::centi::num == 1 && boost::centi::den == 100ULL);
STATIC_ASSERT(boost::deci::num == 1 && boost::deci::den == 10ULL);
STATIC_ASSERT(boost::deca::num == 10ULL && boost::deca::den == 1);
STATIC_ASSERT(boost::hecto::num == 100ULL && boost::hecto::den == 1);
STATIC_ASSERT(boost::kilo::num == 1000ULL && boost::kilo::den == 1);
STATIC_ASSERT(boost::mega::num == 1000000ULL && boost::mega::den == 1);
STATIC_ASSERT(boost::giga::num == 1000000000ULL && boost::giga::den == 1);
STATIC_ASSERT(boost::tera::num == 1000000000000ULL && boost::tera::den == 1);
STATIC_ASSERT(boost::peta::num == 1000000000000000ULL && boost::peta::den == 1);
STATIC_ASSERT(boost::exa::num == 1000000000000000000ULL && boost::exa::den == 1);
