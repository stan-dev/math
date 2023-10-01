// Copyright 2022 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/fraction.hpp>
#include <boost/histogram/utility/jeffreys_interval.hpp>
#include <limits>
#include "is_close.hpp"
#include "throw_exception.hpp"

using namespace boost::histogram::utility;
using namespace boost::histogram::accumulators;

template <class T>
void test() {
  // reference: table A.1 in
  // L.D. Brown, T.T. Cai, A. DasGupta, Statistical Science 16 (2001) 101–133,
  // doi:10.1214/ss/1009213286

  const T atol = 0.001;

  jeffreys_interval<T> iv(confidence_level{0.95});

  {
    auto p = iv(0, 7);
    BOOST_TEST_IS_CLOSE(p.first, 0, atol);
    BOOST_TEST_IS_CLOSE(p.second, 0.41, atol);
  }

  {
    auto p = iv(1, 6);
    BOOST_TEST_IS_CLOSE(p.first, 0, atol);
    BOOST_TEST_IS_CLOSE(p.second, 0.501, atol);
  }

  {
    auto p = iv(2, 5);
    BOOST_TEST_IS_CLOSE(p.first, 0.065, atol);
    BOOST_TEST_IS_CLOSE(p.second, 0.648, atol);
  }

  {
    auto p = iv(3, 4);
    BOOST_TEST_IS_CLOSE(p.first, 0.139, atol);
    BOOST_TEST_IS_CLOSE(p.second, 0.766, atol);
  }

  {
    auto p = iv(4, 7 - 4);
    BOOST_TEST_IS_CLOSE(p.first, 0.234, atol);
    BOOST_TEST_IS_CLOSE(p.second, 0.861, atol);
  }

  // extrapolated from table
  {
    auto p = iv(5, 2);
    BOOST_TEST_IS_CLOSE(p.first, 1 - 0.648, atol);
    BOOST_TEST_IS_CLOSE(p.second, 1 - 0.065, atol);
  }

  // extrapolated from table
  {
    auto p = iv(6, 1);
    BOOST_TEST_IS_CLOSE(p.first, 1 - 0.501, atol);
    BOOST_TEST_IS_CLOSE(p.second, 1, atol);
  }

  // extrapolated from table
  {
    auto p = iv(7, 0);
    BOOST_TEST_IS_CLOSE(p.first, 1 - 0.41, atol);
    BOOST_TEST_IS_CLOSE(p.second, 1, atol);
  }

  // not in table
  {
    auto p = iv(0, 1);
    BOOST_TEST_IS_CLOSE(p.first, 0, atol);
    BOOST_TEST_IS_CLOSE(p.second, 0.975, atol);

    fraction<T> f(0, 1);
    const auto y = iv(f);
    BOOST_TEST_IS_CLOSE(y.first, 0.0, atol);
    BOOST_TEST_IS_CLOSE(y.second, 0.975, atol);
  }

  // not in table
  {
    auto p = iv(1, 0);
    BOOST_TEST_IS_CLOSE(p.first, 0.025, atol);
    BOOST_TEST_IS_CLOSE(p.second, 1, atol);

    fraction<T> f(1, 0);
    const auto y = iv(f);
    BOOST_TEST_IS_CLOSE(y.first, 0.025, atol);
    BOOST_TEST_IS_CLOSE(y.second, 1, atol);
  }
}

int main() {

  test<float>();
  test<double>();

  return boost::report_errors();
}
