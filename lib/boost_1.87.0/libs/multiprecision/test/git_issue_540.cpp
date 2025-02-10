///////////////////////////////////////////////////////////////////////////////
//  Copyright 2023 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include "test.hpp"

int main()
{
   using Integer  = boost::multiprecision::cpp_int;
   using Rational = boost::multiprecision::cpp_rational;

   Integer  one                   = 1;
   Integer  tenThousand           = 10000;
   Integer  threeFourFiveSix      = 3456;
   Rational oneInTenThousand      = Rational(one) / tenThousand;
   Rational oneInThreeFourFiveSix = Rational(one) / threeFourFiveSix;

   Rational result_1(1, 10000);
   Rational result_2(1, 3456);

   BOOST_CHECK_EQUAL(oneInTenThousand, result_1);
   BOOST_CHECK_EQUAL(oneInThreeFourFiveSix, result_2);
}
