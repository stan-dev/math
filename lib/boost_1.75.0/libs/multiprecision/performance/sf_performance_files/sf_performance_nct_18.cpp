///////////////////////////////////////////////////////////////
//  Copyright 2011 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#include "../sf_performance.hpp"

void nct_tests_18()
{
#ifdef TEST_CPP_BIN_FLOAT
   time_proc("Non-central T Distribution (100 digit precision)", "cpp_bin_float_100", test_nct<cpp_bin_float_100>);
#endif
}
