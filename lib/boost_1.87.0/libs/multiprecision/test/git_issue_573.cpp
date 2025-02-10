///////////////////////////////////////////////////////////////////////////////
//  Copyright 2023 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// In C++23 and later std::numeric_limits<>::has_denorm along with type
// float_denorm_style are deprecated, and MSVC generates either a C4996
// warning, or when building with /cdl (the recomended default for the IDE)
// a hard error.  It is sufficient to include our headers to reproduce the
// issue, unless all uses of float_denorm_style are guarded by #pragma warning(disable : 4996)
// 
// Note that other compilers may follow suite in due course, and in time float_denorm_style
// will be removed completely.  For now though we DO need to provide numeric_limits<>::has_denorm
// in our specializations, as it's only deprecated, not removed as of C++23.
//
// See https://github.com/boostorg/multiprecision/issues/573

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>

int main()
{
   boost::multiprecision::cpp_dec_float_50 f1(2);
   std::cout << tgamma(f1) << std::endl;
   boost::multiprecision::cpp_bin_float_50 f2(2);
   std::cout << tgamma(f2) << std::endl;
   boost::multiprecision::cpp_int i = 2;
   std::cout << i * i << std::endl;
}

