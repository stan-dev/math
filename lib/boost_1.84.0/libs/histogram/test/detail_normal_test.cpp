// Copyright 2022 Hans Dembinski, Jay Gohil
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/detail/normal.hpp>
#include <boost/math/distributions/normal.hpp>
#include <limits>
#include "is_close.hpp"
#include "throw_exception.hpp"

using namespace boost::histogram::detail;
namespace bm = boost::math;

int main() {
  const auto x = {-6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 6.0};

  const double eps = std::numeric_limits<double>::epsilon();
  bm::normal norm;

  // cdf
  {
    for (auto&& xi : x) {
      const double expected = bm::cdf(norm, xi);
      BOOST_TEST_IS_CLOSE(normal_cdf(xi), expected, eps);
    }
  }

  // ppf
  {
    for (auto&& xi : x) {
      const double p = bm::cdf(norm, xi);
      BOOST_TEST_IS_CLOSE(normal_ppf(p), xi, 1e-8);
    }
  }

  return boost::report_errors();
}
