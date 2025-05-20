///////////////////////////////////////////////////////////////
//  Copyright 2013 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

// Demonstrations of using Boost.Multiprecision float128 quad type.
// (Only available using GCC compiler).

// Contains Quickbook markup in comments.

//[float128_eg
#include <boost/type_index.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>

int main()
{
   using namespace boost::multiprecision; // Potential to cause name collisions?
   // using boost::multiprecision::float128; // is safer.

/*`The type float128 provides operations at 128-bit precision with
[@https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format#IEEE_754_quadruple-precision_binary_floating-point_format:_binary128 Quadruple-precision floating-point format]
and have full `std::numeric_limits` support:
*/
   float128 b = 2;
//` There are  15 bits of (biased) binary exponent and 113-bits of significand precision
   std::cout << std::numeric_limits<float128>::digits << std::endl;
//` or 33 decimal places:
   std::cout << std::numeric_limits<float128>::digits10 << std::endl;
   //` We can use any C++ std library function, so let's show all the at-most 36 potentially significant digits, and any trailing zeros, as well:
    std::cout.setf(std::ios_base::showpoint); // Include any trailing zeros.
   std::cout << std::setprecision(std::numeric_limits<float128>::max_digits10)
      << log(b) << std::endl; // Shows log(2) = 0.693147180559945309417232121458176575
   //` We can also use any function from Boost.Math, for example, the 'true gamma' function `tgamma`:
   std::cout << boost::math::tgamma(b) << std::endl;
   /*` And since we have an extended exponent range, we can generate some really large
    numbers here (4.02387260077093773543702433923004111e+2564):
   */
   std::cout << boost::math::tgamma(float128(1000)) << std::endl;
   /*` We can declare constants using GCC or Intel's native types, and literals with the Q suffix, and these can be declared `constexpr` if required:
   */
   // std::numeric_limits<float128>::max_digits10 = 36
   constexpr float128 pi = 3.14159265358979323846264338327950288Q;
   std::cout.precision(std::numeric_limits<float128>::max_digits10);
   std::cout << "pi = " << pi << std::endl; //pi = 3.14159265358979323846264338327950280
//] [/float128_eg]

   // Full generation for numeric limits below:
   using boost::multiprecision::float128;
   std::cout << std::boolalpha;
   std::cout << std::setprecision(std::numeric_limits<float128>::max_digits10);

   std::cout << "GCC " << __VERSION__ << std::endl << std::endl;
   std::cout << "Type name: boost::multiprecision::float128" << std::endl;
   std::cout << "Full name: " << boost::typeindex::type_id<float128>().pretty_name() << std::endl << std::endl;
#if defined(__cpp_lib_type_trait_variable_templates) && (__cpp_lib_type_trait_variable_templates >= 201510L)
   std::cout << "std::is_fundamental<> = " << std::is_fundamental_v<float128> << std::endl;
   std::cout << "boost::multiprecision::detail::is_signed<> = " << boost::multiprecision::detail::is_signed_v<float128> << std::endl;
   std::cout << "boost::multiprecision::detail::is_unsigned<> = " << boost::multiprecision::detail::is_unsigned_v<float128> << std::endl;
   std::cout << "boost::multiprecision::detail::is_integral<> = " << boost::multiprecision::detail::is_integral_v<float128> << std::endl;
   std::cout << "boost::multiprecision::detail::is_arithmetic<> = " << boost::multiprecision::detail::is_arithmetic_v<float128> << std::endl;
   std::cout << "std::is_const<> = " << std::is_const_v<float128> << std::endl;
   std::cout << "std::is_trivial<> = " << std::is_trivial_v<float128> << std::endl;
   std::cout << "std::is_trivially_copyable<> = " << std::is_trivially_copyable_v<float128> << std::endl;
   std::cout << "std::is_standard_layout<> = " << std::is_standard_layout_v<float128> << std::endl;
   std::cout << "std::is_pod<> = " << std::is_pod_v<float128> << std::endl;
#endif
   std::cout << "std::numeric_limits<>::is_exact = " << std::numeric_limits<float128>::is_exact << std::endl;
   std::cout << "std::numeric_limits<>::is_bounded = " << std::numeric_limits<float128>::is_bounded << std::endl;
   std::cout << "std::numeric_limits<>::is_modulo = " << std::numeric_limits<float128>::is_modulo << std::endl;
   std::cout << "std::numeric_limits<>::is_iec559 = " << std::numeric_limits<float128>::is_iec559 << std::endl;
   std::cout << "std::numeric_limits<>::traps = " << std::numeric_limits<float128>::traps << std::endl;
   std::cout << "std::numeric_limits<>::tinyness_before = " << std::numeric_limits<float128>::tinyness_before << std::endl;
   std::cout << "std::numeric_limits<>::max() = " << (std::numeric_limits<float128>::max)() << std::endl;
   std::cout << "std::numeric_limits<>::min() = " << (std::numeric_limits<float128>::min)() << std::endl;
   std::cout << "std::numeric_limits<>::lowest() = " << std::numeric_limits<float128>::lowest() << std::endl;
   std::cout << "std::numeric_limits<>::min_exponent = " << std::numeric_limits<float128>::min_exponent << std::endl;
   std::cout << "std::numeric_limits<>::max_exponent = " << std::numeric_limits<float128>::max_exponent << std::endl;
   std::cout << "std::numeric_limits<>::epsilon() = " << std::numeric_limits<float128>::epsilon() << std::endl;
   std::cout << "std::numeric_limits<>::radix = " << std::numeric_limits<float128>::radix << std::endl;
   std::cout << "std::numeric_limits<>::digits = " << std::numeric_limits<float128>::digits << std::endl;
   std::cout << "std::numeric_limits<>::max_digits10 = " << std::numeric_limits<float128>::max_digits10 << std::endl;
   std::cout << "std::numeric_limits<>::has_denorm = " << std::numeric_limits<float128>::has_denorm << std::endl;
   std::cout << "std::numeric_limits<>::denorm_min() = " << std::numeric_limits<float128>::denorm_min() << std::endl;
   std::cout << "std::numeric_limits<>::has_denorm_loss = " << std::numeric_limits<float128>::has_denorm_loss << std::endl;
   std::cout << "std::numeric_limits<>::has_signaling_NaN = " << std::numeric_limits<float128>::has_signaling_NaN << std::endl;
   std::cout << "std::numeric_limits<>::quiet_NaN() = " << std::numeric_limits<float128>::quiet_NaN() << std::endl;
   std::cout << "std::numeric_limits<>::infinity() = " << std::numeric_limits<float128>::infinity() << std::endl;

   return 0;
}

/*
//[float128_numeric_limits

GCC 14.1.0

Type name: boost::multiprecision::float128
Full name: boost::multiprecision::number<boost::multiprecision::backends::float128_backend, (boost::multiprecision::expression_template_option)0>

std::is_fundamental<> = false
boost::multiprecision::detail::is_signed<> = false
boost::multiprecision::detail::is_unsigned<> = false
boost::multiprecision::detail::is_integral<> = false
boost::multiprecision::detail::is_arithmetic<> = false
std::is_const<> = false
std::is_trivial<> = false
std::is_trivially_copyable<> = true
std::is_standard_layout<> = true
std::is_pod<> = false
std::numeric_limits<>::is_exact = false
std::numeric_limits<>::is_bounded = true
std::numeric_limits<>::is_modulo = false
std::numeric_limits<>::is_iec559 = true
std::numeric_limits<>::traps = false
std::numeric_limits<>::tinyness_before = false
std::numeric_limits<>::max() = 1.18973149535723176508575932662800702e+4932
std::numeric_limits<>::min() = 3.3621031431120935062626778173217526e-4932
std::numeric_limits<>::lowest() = -1.18973149535723176508575932662800702e+4932
std::numeric_limits<>::min_exponent = -16381
std::numeric_limits<>::max_exponent = 16384
std::numeric_limits<>::epsilon() = 1.92592994438723585305597794258492732e-34
std::numeric_limits<>::radix = 2
std::numeric_limits<>::digits = 113
std::numeric_limits<>::max_digits10 = 36
std::numeric_limits<>::has_denorm = 1
std::numeric_limits<>::denorm_min() = 6.47517511943802511092443895822764655e-4966
std::numeric_limits<>::has_denorm_loss = true
std::numeric_limits<>::has_signaling_NaN = false
std::numeric_limits<>::quiet_NaN() = nan
std::numeric_limits<>::infinity() = inf

//] [/float128_numeric_limits]
*/


