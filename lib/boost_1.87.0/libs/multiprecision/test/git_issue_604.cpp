///////////////////////////////////////////////////////////////////////////////
//  Copyright 2024 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#if !defined(TEST_MPF_50) && !defined(TEST_CPP_DEC_FLOAT) && !defined(TEST_MPFR_50) && !defined(TEST_FLOAT128) && !defined(TEST_CPP_BIN_FLOAT)
#define TEST_MPF_50
#define TEST_MPFR_50
#define TEST_CPP_DEC_FLOAT
#define TEST_FLOAT128
#define TEST_CPP_BIN_FLOAT

#ifdef _MSC_VER
#pragma message("CAUTION!!: No backend type specified so testing everything.... this will take some time!!")
#endif
#ifdef __GNUC__
#pragma warning "CAUTION!!: No backend type specified so testing everything.... this will take some time!!"
#endif

#endif

#if defined(TEST_MPF_50)
#include <boost/multiprecision/gmp.hpp>
#endif
#if defined(TEST_MPFR_50)
#include <boost/multiprecision/mpfr.hpp>
#endif
#ifdef TEST_CPP_DEC_FLOAT
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif
#ifdef TEST_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif
#ifdef TEST_CPP_BIN_FLOAT
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

#include <iostream>
#include "test.hpp"

template <class T>
void test()
{
   {
      const T val1{6};
      const T val2{3};
      T       r = boost::multiprecision::fmod(val1, val2);
      BOOST_CHECK_EQUAL(r, T(0));
   }
   {
      const T val1{-6};
      const T val2{3};
      T       r = boost::multiprecision::fmod(val1, val2);
      BOOST_CHECK_EQUAL(r, T(0));
   }
   {
      const T val1{6};
      const T val2{-3};
      T       r = boost::multiprecision::fmod(val1, val2);
      BOOST_CHECK_EQUAL(r, T(0));
   }
   {
      const T val1{-6};
      const T val2{-3};
      T       r = boost::multiprecision::fmod(val1, val2);
      BOOST_CHECK_EQUAL(r, T(0));
   }
}



int main()
{
   using namespace boost::multiprecision;

#ifdef TEST_MPF_50
   test<boost::multiprecision::mpf_float_50>();
   test<boost::multiprecision::mpf_float_100>();
#endif
#ifdef TEST_MPFR_50
   test<boost::multiprecision::mpfr_float_50>();
   test<boost::multiprecision::mpfr_float_100>();
#endif
#ifdef TEST_CPP_DEC_FLOAT
   test<boost::multiprecision::cpp_dec_float_50>();
   test<boost::multiprecision::cpp_dec_float_100>();
#endif
#ifdef TEST_FLOAT128
   test<boost::multiprecision::float128>();
#endif
#ifdef TEST_CPP_BIN_FLOAT
   test<boost::multiprecision::cpp_bin_float_50>();
   test<boost::multiprecision::number<boost::multiprecision::cpp_bin_float<35, boost::multiprecision::digit_base_10, std::allocator<char>, long long> > >();
#endif
   return boost::report_errors();
}
