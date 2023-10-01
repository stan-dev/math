// Copyright 2018 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_HISTOGRAM_TEST_IS_CLOSE_HPP
#define BOOST_HISTOGRAM_TEST_IS_CLOSE_HPP

#include <boost/core/lightweight_test.hpp>
#include <cmath>

struct is_close {
  bool operator()(double a, double b) { return std::abs(a - b) < atol; }
  double atol;
};

#define BOOST_TEST_IS_CLOSE(a, b, atol) BOOST_TEST_WITH(a, b, is_close{atol})

#endif
