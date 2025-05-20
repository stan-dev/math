///////////////////////////////////////////////////////////////
//  Copyright 2012 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#include <boost/multiprecision/cpp_bin_float.hpp>

#include "libs/multiprecision/test/test_arithmetic.hpp"

#if 0
template <unsigned Digits, boost::multiprecision::backends::digit_base_type DigitBase, class Allocator, class Exponent, Exponent MinExponent, Exponent MaxExponent, boost::multiprecision::expression_template_option ET>
struct related_type<boost::multiprecision::number<boost::multiprecision::cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinExponent, MaxExponent>, ET> >
{
   typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinExponent, MaxExponent>, ET>                                                                                                            number_type;
   typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<((std::numeric_limits<number_type>::digits / 2) > std::numeric_limits<long double>::digits ? Digits / 2 : Digits), DigitBase, Allocator, Exponent, MinExponent, MaxExponent>, ET> type;
};
#endif

int main()
{
   using namespace boost::multiprecision;
   typedef number<backends::cpp_bin_float<11, backends::digit_base_2, void, std::int16_t, -14, 15>, et_off> float16_t;
   typedef number<backends::cpp_bin_float<8, backends::digit_base_2, void, std::int16_t, -126, 127>, et_off> bfloat16_t;
   
   test<float16_t>();
   test<bfloat16_t>();
   return boost::report_errors();
}
