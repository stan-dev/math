// Copyright 2015-2018 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/ostream.hpp>
#include <boost/histogram/accumulators/weighted_mean.hpp>
#include <boost/histogram/weight.hpp>
#include <sstream>
#include "is_close.hpp"
#include "throw_exception.hpp"
#include "utility_str.hpp"

using namespace boost::histogram;
using namespace std::literals;

int main() {
  using m_t = accumulators::weighted_mean<double>;
  m_t a;
  BOOST_TEST_EQ(a.sum_of_weights(), 0);
  BOOST_TEST_EQ(a, m_t{});

  a(weight(0.5), 1);
  a(weight(1.0), 2);
  a(weight(0.5), 3);

  BOOST_TEST_EQ(a.sum_of_weights(), 2);
  BOOST_TEST_EQ(a.sum_of_weights_squared(), 1.5);
  BOOST_TEST_EQ(a.value(), 2);
  BOOST_TEST_IS_CLOSE(a.variance(), 0.8, 1e-3);

  BOOST_TEST_EQ(str(a), "weighted_mean(2, 2, 0.8)"s);
  BOOST_TEST_EQ(str(a, 25, false), " weighted_mean(2, 2, 0.8)"s);
  BOOST_TEST_EQ(str(a, 25, true), "weighted_mean(2, 2, 0.8) "s);

  auto b = a;
  b += a; // same as feeding all samples twice

  BOOST_TEST_EQ(b.sum_of_weights(), 4);
  BOOST_TEST_EQ(b.value(), 2);
  BOOST_TEST_IS_CLOSE(b.variance(), 0.615, 1e-3);

  BOOST_TEST_EQ(m_t() += m_t(), m_t());
  BOOST_TEST_EQ(m_t(1, 2, 3, 4) += m_t(), m_t(1, 2, 3, 4));
  BOOST_TEST_EQ(m_t() += m_t(1, 2, 3, 4), m_t(1, 2, 3, 4));

  return boost::report_errors();
}
