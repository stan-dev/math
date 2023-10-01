// Copyright 2019 Henry Schreiner, Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// The header windows.h and possibly others illegially use unprotected defines
#define small macro_substitution_this_shouldn_t_be_happening
#define min(A, B) macro_substitution_this_shouldn_t_be_happening
#define max(A, B) macro_substitution_this_shouldn_t_be_happening
// which violates the C++ standard. We make sure here that including our headers work
// nevertheless by avoiding these preprocessing tokens or by preventing their macro
// substitution. For more details, see https://github.com/boostorg/histogram/issues/342

// include all Boost.Histogram header here; see odr_main_test.cpp for details
#include <boost/histogram.hpp>
#include <boost/histogram/ostream.hpp>
#include <boost/histogram/serialization.hpp>

#include <boost/histogram/detail/ignore_deprecation_warning_begin.hpp>
#include <boost/histogram/detail/ignore_deprecation_warning_end.hpp>
#include <boost/histogram/utility/clopper_pearson_interval.hpp>
#include <boost/histogram/utility/jeffreys_interval.hpp>
#include <boost/histogram/utility/wald_interval.hpp>
