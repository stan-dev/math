///////////////////////////////////////////////////////////////////////////////
//  Copyright 2020 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/static_assert.hpp>
#include "test.hpp"

using namespace boost::multiprecision;

int main()
{
   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_single>::digits10, 6);
   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_single>::max_digits10, 9);

   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_double>::digits10, 15);
   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_double>::max_digits10, 17);

   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_double_extended>::digits10, 18);
   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_double_extended>::max_digits10, 21);

   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_quad>::digits10, 33);
   BOOST_CHECK_EQUAL(std::numeric_limits<cpp_bin_float_quad>::max_digits10, 36);

   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_single>::digits10 == 6);
   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_single>::max_digits10 == 9);

   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_double>::digits10 == 15);
   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_double>::max_digits10 == 17);

   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_double_extended>::digits10 == 18);
   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_double_extended>::max_digits10 == 21);

   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_quad>::digits10 == 33);
   BOOST_STATIC_ASSERT(std::numeric_limits<cpp_bin_float_quad>::max_digits10 == 36);

   return boost::report_errors();
}
