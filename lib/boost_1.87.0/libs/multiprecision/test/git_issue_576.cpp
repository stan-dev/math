///////////////////////////////////////////////////////////////////////////////
//  Copyright 2023 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <iostream>
#include "test.hpp"

template <class Float, class F>
void test(F f, bool big)
{
   Float f2(f);
   if (big)
   {
      F tol = f * static_cast<F>(std::numeric_limits<Float>::epsilon());
      F max  = static_cast<F>((std::numeric_limits<Float>::max)());
      F diff = static_cast<F>(f2) - f;
      if (diff < 0)
         diff = -diff;
      if (f > max)
      {
         BOOST_CHECK(isinf(f2) || (diff < tol));
      }
      else
      {
         BOOST_CHECK_LE(diff, tol);
      }
   }
   else
   {
      BOOST_CHECK_EQUAL(static_cast<F>(f2), f);
   }
}


template <class Float>
void test()
{
   std::int64_t i = static_cast<std::int64_t>(1) << (std::numeric_limits<Float>::max_exponent + 1);
   Float        f(i);
   BOOST_CHECK_EQUAL(f, std::numeric_limits<Float>::infinity());
   
   if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<float>::digits)
   {
      BOOST_CHECK_EQUAL(Float(static_cast<float>(i)), std::numeric_limits<Float>::infinity());
   }
   if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<double>::digits)
   {
      BOOST_CHECK_EQUAL(Float(static_cast<double>(i)), std::numeric_limits<Float>::infinity());
   }
   if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<long double>::digits)
   {
      BOOST_CHECK_EQUAL(Float(static_cast<long double>(i)), std::numeric_limits<Float>::infinity());
   }
#ifdef BOOST_HAS_FLOAT128
   if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<__float128>::digits)
   {
      BOOST_CHECK_EQUAL(Float(static_cast<__float128>(i)), std::numeric_limits<Float>::infinity());
   }
#endif
   
   --i;

   while (i)
   {
      Float f2(i);
      BOOST_CHECK_NE(f2, std::numeric_limits<Float>::infinity());
      bool big = boost::multiprecision::msb(i) >= std::numeric_limits<Float>::digits;
      if (big)
      {
         BOOST_CHECK_LE(static_cast<std::int64_t>(f2), i);
      }
      else
      {
         BOOST_CHECK_EQUAL(static_cast<std::int64_t>(f2), i);
      }

      if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<float>::digits)
      {
         test<Float>(static_cast<float>(i), big);
         for (int exp = -1; exp >= std::numeric_limits<Float>::min_exponent; --exp)
            test<Float>(std::ldexp(static_cast<float>(i), exp), big);
      }
      if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<double>::digits)
      {
         test<Float>(static_cast<double>(i), big);
         for (int exp = -1; exp >= std::numeric_limits<Float>::min_exponent; --exp)
            test<Float>(std::ldexp(static_cast<double>(i), exp), big);
      }
      if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<long double>::digits)
      {
         test<Float>(static_cast<long double>(i), big);
         for (int exp = -1; exp >= std::numeric_limits<Float>::min_exponent; --exp)
            test<Float>(std::ldexp(static_cast<long double>(i), exp), big);
      }
#ifdef BOOST_HAS_FLOAT128
      if (std::numeric_limits<Float>::max_exponent < std::numeric_limits<__float128>::digits)
      {
         test<Float>(static_cast<__float128>(i), big);
         for (int exp = -1; exp >= std::numeric_limits<Float>::min_exponent; --exp)
            test<Float>(ldexpq(static_cast<__float128>(i), exp), big);
      }
#endif

      --i;
   }
}



int main()
{
   using namespace boost::multiprecision;
   typedef number<backends::cpp_bin_float<11, backends::digit_base_2, void, std::int16_t, -14, 15>, et_off>  float16_t;
   //typedef number<backends::cpp_bin_float<8, backends::digit_base_2, void, std::int16_t, -126, 127>, et_off> bfloat16_t;

   test<float16_t>();
   //test<bfloat16_t>();

   return boost::report_errors();
}
