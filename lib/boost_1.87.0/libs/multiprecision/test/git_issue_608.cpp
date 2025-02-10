///////////////////////////////////////////////////////////////////////////////
//  Copyright 2024 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <functional>
#include <boost/multiprecision/cpp_bin_float.hpp>

using big_float_type = boost::multiprecision::cpp_bin_float_100;

int main()
{
   static_assert(boost::multiprecision::is_compatible_arithmetic_type<std::reference_wrapper<big_float_type>, big_float_type>::value == 0, "This should not be a compatible type");
}

