// Copyright 2022 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/utility/binomial_proportion_interval.hpp>
#include "is_close.hpp"
#include "throw_exception.hpp"

using namespace boost::histogram::utility;

int main() {

  // confidence_level <-> deviation
  {
    confidence_level cl1 = deviation{1};
    BOOST_TEST_IS_CLOSE(static_cast<double>(cl1), 0.683, 1e-3);
    BOOST_TEST_IS_CLOSE(static_cast<double>(deviation(cl1)), 1, 1e-8);
    confidence_level cl2 = deviation{2};
    BOOST_TEST_IS_CLOSE(static_cast<double>(cl2), 0.954, 1e-3);
    BOOST_TEST_IS_CLOSE(static_cast<double>(deviation(cl2)), 2, 1e-8);
    confidence_level cl3 = deviation{3};
    BOOST_TEST_IS_CLOSE(static_cast<double>(cl3), 0.997, 1e-3);
    BOOST_TEST_IS_CLOSE(static_cast<double>(deviation(cl3)), 3, 1e-8);
  }

  return boost::report_errors();
}
