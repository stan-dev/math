///////////////////////////////////////////////////////////////////////////////
//  Copyright 2023 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/multiprecision/cpp_complex.hpp>

#include "eigen.hpp"

int main()
{
   using namespace boost::multiprecision;
   test_complex_type<cpp_complex_50>();
   return 0;
}
