///////////////////////////////////////////////////////////////////////////////
//  Copyright 2023 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/complex_adaptor.hpp>
#include <boost/multiprecision/logged_adaptor.hpp>
#include <boost/multiprecision/debug_adaptor.hpp>
#include <iostream>
#include "test.hpp"

template <class T>
void test()
{
   T val;
   unsigned n = T::default_precision();
   T::default_precision(n);
   n = T::thread_default_precision();
   T::thread_default_precision(n);
   n = val.precision();
   val.precision(n);
}



int main()
{
   using namespace boost::multiprecision;

   test<mpfr_float>();

   using c = number<complex_adaptor<mpfr_float::backend_type>, et_on>;

   test<c>();

   using l = number<logged_adaptor<mpfr_float::backend_type>, et_on>;

   test<l>();

   using d = number<debug_adaptor<mpfr_float::backend_type>, et_on>;

   test<d>();

   return boost::report_errors();
}
