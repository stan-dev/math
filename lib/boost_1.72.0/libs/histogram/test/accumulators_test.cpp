// Copyright 2015-2018 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/mean.hpp>
#include <boost/histogram/accumulators/ostream.hpp>
#include <boost/histogram/accumulators/sum.hpp>
#include <boost/histogram/accumulators/thread_safe.hpp>
#include <boost/histogram/accumulators/weighted_mean.hpp>
#include <boost/histogram/accumulators/weighted_sum.hpp>
#include <boost/histogram/weight.hpp>
#include <sstream>
#include "is_close.hpp"
#include "throw_exception.hpp"

using namespace boost::histogram;
using namespace std::literals;

template <class T>
auto str(const T& t, int w = 0, bool left = true) {
  std::ostringstream os;
  os.width(w);
  if (left)
    os << std::left;
  else
    os << std::right;
  os << t;
  return os.str();
}

int main() {
  {
    using w_t = accumulators::weighted_sum<double>;
    w_t w;
    BOOST_TEST_EQ(str(w), "weighted_sum(0, 0)"s);
    BOOST_TEST_EQ(str(w, 20, false), "  weighted_sum(0, 0)"s);
    BOOST_TEST_EQ(str(w, 20, true), "weighted_sum(0, 0)  "s);
    BOOST_TEST_EQ(w, w_t{});

    BOOST_TEST_EQ(w, w_t(0));
    BOOST_TEST_NE(w, w_t(1));
    w = w_t(1);
    BOOST_TEST_EQ(w.value(), 1);
    BOOST_TEST_EQ(w.variance(), 1);
    BOOST_TEST_EQ(w, 1);
    BOOST_TEST_NE(w, 2);

    w += 2;
    BOOST_TEST_EQ(w.value(), 3);
    BOOST_TEST_EQ(w.variance(), 5);
    BOOST_TEST_EQ(w, w_t(3, 5));
    BOOST_TEST_NE(w, w_t(3));

    w += w_t(1, 2);
    BOOST_TEST_EQ(w.value(), 4);
    BOOST_TEST_EQ(w.variance(), 7);

    // consistency: a weighted counter increased by weight 1 multiplied
    // by 2 must be the same as a weighted counter increased by weight 2
    w_t u(0);
    ++u;
    u *= 2;
    BOOST_TEST_EQ(u, w_t(2, 4));

    w_t v(0);
    v += 2;
    BOOST_TEST_EQ(u, v);

    // conversion to RealType
    w_t y(1, 2);
    BOOST_TEST_NE(y, 1);
    BOOST_TEST_EQ(static_cast<double>(y), 1);

    BOOST_TEST_EQ(w_t() += w_t(), w_t());
  }

  {
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
  }

  {
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
  }

  {
    double bad_sum = 0;
    bad_sum += 1;
    bad_sum += 1e100;
    bad_sum += 1;
    bad_sum += -1e100;
    BOOST_TEST_EQ(bad_sum, 0); // instead of 2

    using s_t = accumulators::sum<double>;
    s_t sum;
    ++sum;
    BOOST_TEST_EQ(sum.large(), 1);
    BOOST_TEST_EQ(sum.small(), 0);
    BOOST_TEST_EQ(str(sum), "sum(1 + 0)"s);
    BOOST_TEST_EQ(str(sum, 15, false), "     sum(1 + 0)"s);
    BOOST_TEST_EQ(str(sum, 15, true), "sum(1 + 0)     "s);

    sum += 1e100;
    BOOST_TEST_EQ(str(sum), "sum(1e+100 + 1)"s);
    ++sum;
    BOOST_TEST_EQ(str(sum), "sum(1e+100 + 2)"s);
    sum += -1e100;
    BOOST_TEST_EQ(str(sum), "sum(0 + 2)"s);
    BOOST_TEST_EQ(sum, 2); // correct answer
    BOOST_TEST_EQ(sum.large(), 0);
    BOOST_TEST_EQ(sum.small(), 2);

    accumulators::sum<double> a(3), b(2), c(3);
    BOOST_TEST_LT(b, c);
    BOOST_TEST_LE(b, c);
    BOOST_TEST_LE(a, c);
    BOOST_TEST_GT(a, b);
    BOOST_TEST_GE(a, b);
    BOOST_TEST_GE(a, c);

    BOOST_TEST_EQ(s_t() += s_t(), s_t());
  }

  {
    using s_t = accumulators::weighted_sum<accumulators::sum<double>>;
    s_t w;

    ++w;
    w += 1e100;
    ++w;
    w += -1e100;

    BOOST_TEST_EQ(w.value(), 2);
    BOOST_TEST_EQ(w.variance(), 2e200);

    BOOST_TEST_EQ(s_t() += s_t(), s_t());
  }

  {
    using ts_t = accumulators::thread_safe<int>;
    ts_t i;
    ++i;
    i += 1000;

    BOOST_TEST_EQ(i, 1001);
    BOOST_TEST_EQ(str(i), "1001"s);

    BOOST_TEST_EQ(ts_t() += ts_t(), ts_t());
  }

  return boost::report_errors();
}
