// Copyright 2015-2018 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/fraction.hpp>
#include <boost/histogram/accumulators/ostream.hpp>
#include <boost/histogram/utility/wilson_interval.hpp>
#include <cmath>
#include <limits>
#include "is_close.hpp"
#include "str.hpp"
#include "throw_exception.hpp"

using namespace boost::histogram;
using namespace std::literals;

template <class T>
void run_tests() {
  using f_t = accumulators::fraction<T>;

  const double eps = std::numeric_limits<typename f_t::real_type>::epsilon();

  {
    f_t f;
    BOOST_TEST_EQ(f.successes(), 0);
    BOOST_TEST_EQ(f.failures(), 0);
    BOOST_TEST(std::isnan(f.value()));
    BOOST_TEST(std::isnan(f.variance()));

    const auto ci = f.confidence_interval();
    BOOST_TEST(std::isnan(ci.first));
    BOOST_TEST(std::isnan(ci.second));
  }

  {
    f_t f;
    f(true);
    BOOST_TEST_EQ(f.successes(), 1);
    BOOST_TEST_EQ(f.failures(), 0);
    BOOST_TEST_EQ(str(f), "fraction(1, 0)"s);
    f(false);
    BOOST_TEST_EQ(f.successes(), 1);
    BOOST_TEST_EQ(f.failures(), 1);
    BOOST_TEST_EQ(str(f), "fraction(1, 1)"s);
    BOOST_TEST_EQ(str(f, 20), "fraction(1, 1)      "s);
  }

  {
    f_t f(3, 1);
    BOOST_TEST_EQ(f.successes(), 3);
    BOOST_TEST_EQ(f.failures(), 1);
    BOOST_TEST_EQ(f.value(), 0.75);
    BOOST_TEST_IS_CLOSE(f.variance(), 0.75 * (1 - 0.75) / 4, eps);

    const auto ci = f.confidence_interval();
    const auto expected = utility::wilson_interval<double>()(3, 1);
    BOOST_TEST_IS_CLOSE(ci.first, expected.first, eps);
    BOOST_TEST_IS_CLOSE(ci.second, expected.second, eps);
  }

  {
    f_t f(0, 1);
    BOOST_TEST_EQ(f.successes(), 0);
    BOOST_TEST_EQ(f.failures(), 1);
    BOOST_TEST_EQ(f.value(), 0);
    BOOST_TEST_EQ(f.variance(), 0);

    const auto ci = f.confidence_interval();
    const auto expected = utility::wilson_interval<double>()(0, 1);
    BOOST_TEST_IS_CLOSE(ci.first, expected.first, eps);
    BOOST_TEST_IS_CLOSE(ci.second, expected.second, eps);
  }

  {
    f_t f(1, 0);
    BOOST_TEST_EQ(f.successes(), 1);
    BOOST_TEST_EQ(f.failures(), 0);
    BOOST_TEST_EQ(f.value(), 1);
    BOOST_TEST_EQ(f.variance(), 0);

    const auto ci = f.confidence_interval();
    const auto expected = utility::wilson_interval<double>()(1, 0);
    BOOST_TEST_IS_CLOSE(ci.first, expected.first, eps);
    BOOST_TEST_IS_CLOSE(ci.second, expected.second, eps);
  }

  {
    f_t a(1, 0), b(0, 1);

    a += b;
    BOOST_TEST_EQ(a, f_t(1, 1));
  }
}

int main() {

  run_tests<int>();
  run_tests<double>();
  run_tests<float>();

  {
    using f_t1 = accumulators::fraction<double>;
    using f_t2 = accumulators::fraction<int>;
    f_t1 f1(5, 3);
    f_t2 f2(f1);
    BOOST_TEST_EQ(f2.successes(), 5);
    BOOST_TEST_EQ(f2.failures(), 3);
  }

  return boost::report_errors();
}
