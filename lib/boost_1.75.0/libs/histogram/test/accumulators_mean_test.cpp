// Copyright 2015-2018 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/mean.hpp>
#include <boost/histogram/accumulators/ostream.hpp>
#include <boost/histogram/weight.hpp>
#include <sstream>
#include "is_close.hpp"
#include "throw_exception.hpp"
#include "utility_str.hpp"

using namespace boost::histogram;
using namespace std::literals;

int main() {
  using m_t = accumulators::mean<double>;
  m_t a;
  BOOST_TEST_EQ(a.count(), 0);
  BOOST_TEST_EQ(a, m_t{});

  a(4);
  a(7);
  a(13);
  a(16);

  BOOST_TEST_EQ(a.count(), 4);
  BOOST_TEST_EQ(a.value(), 10);
  BOOST_TEST_EQ(a.variance(), 30);

  BOOST_TEST_EQ(str(a), "mean(4, 10, 30)"s);
  BOOST_TEST_EQ(str(a, 20, false), "     mean(4, 10, 30)"s);
  BOOST_TEST_EQ(str(a, 20, true), "mean(4, 10, 30)     "s);

  m_t b;
  b(1e8 + 4);
  b(1e8 + 7);
  b(1e8 + 13);
  b(1e8 + 16);

  BOOST_TEST_EQ(b.count(), 4);
  BOOST_TEST_EQ(b.value(), 1e8 + 10);
  BOOST_TEST_EQ(b.variance(), 30);

  auto c = a;
  c += a; // same as feeding all samples twice

  BOOST_TEST_EQ(c.count(), 8);
  BOOST_TEST_EQ(c.value(), 10);
  BOOST_TEST_IS_CLOSE(c.variance(), 25.714, 1e-3);

  // also same as feeding all samples twice
  m_t d;
  d(weight(2), 4);
  d(weight(2), 7);
  d(weight(2), 13);
  d(weight(2), 16);

  BOOST_TEST_EQ(d, c);

  BOOST_TEST_EQ(m_t() += m_t(), m_t());
  BOOST_TEST_EQ(m_t(1, 2, 3) += m_t(), m_t(1, 2, 3));
  BOOST_TEST_EQ(m_t() += m_t(1, 2, 3), m_t(1, 2, 3));

  return boost::report_errors();
}
