///////////////////////////////////////////////////////////////////////////////
//  Copyright 2023 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MP_STANDALONE

#include <iostream>
#include <iomanip>
#include <boost/multiprecision/cpp_bin_float.hpp>

int main()
{
   uint64_t                                  doubleAsInt = 0xBFDFFFFFFFFFFFFC;
   const auto                                value       = std::bit_cast<double>(doubleAsInt);
   boost::multiprecision::cpp_bin_float_quad v           = value;
   return 0;
}
