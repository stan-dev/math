// Copyright 2022 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/fraction.hpp>
#include <boost/histogram/utility/wilson_interval.hpp>
#include <cmath>
#include <limits>
#include "is_close.hpp"
#include "throw_exception.hpp"

using namespace boost::histogram::utility;
using namespace boost::histogram::accumulators;

template <class T>
void test() {
  const double atol = std::numeric_limits<T>::epsilon();

  wilson_interval<T> iv;

  {
    const auto x = iv(0, 0);
    BOOST_TEST(std::isnan(x.first));
    BOOST_TEST(std::isnan(x.second));
  }

  {
    const auto x = iv(0, 1);
    BOOST_TEST_IS_CLOSE(x.first, 0.0, atol);
    BOOST_TEST_IS_CLOSE(x.second, 0.5, atol);

    fraction<T> f(0, 1);
    const auto y = iv(f);
    BOOST_TEST_IS_CLOSE(y.first, 0.0, atol);
    BOOST_TEST_IS_CLOSE(y.second, 0.5, atol);
  }

  {
    const auto x = iv(1, 0);
    BOOST_TEST_IS_CLOSE(x.first, 0.5, atol);
    BOOST_TEST_IS_CLOSE(x.second, 1.0, atol);

    fraction<T> f(1, 0);
    const auto y = iv(f);
    BOOST_TEST_IS_CLOSE(y.first, 0.5, atol);
    BOOST_TEST_IS_CLOSE(y.second, 1.0, atol);
  }

  {
    const auto x = iv(5, 5);
    BOOST_TEST_IS_CLOSE(x.first, 0.3492443277111182, atol);
    BOOST_TEST_IS_CLOSE(x.second, 0.6507556722888819, atol);
  }

  {
    const auto x = iv(1, 9);
    BOOST_TEST_IS_CLOSE(x.first, 0.03887449732033081, atol);
    BOOST_TEST_IS_CLOSE(x.second, 0.23385277540694188, atol);
  }

  {
    const auto x = iv(9, 1);
    BOOST_TEST_IS_CLOSE(x.first, 0.7661472245930581, atol);
    BOOST_TEST_IS_CLOSE(x.second, 0.9611255026796692, atol);
  }
}

int main() {

  test<float>();
  test<double>();

  return boost::report_errors();
}
