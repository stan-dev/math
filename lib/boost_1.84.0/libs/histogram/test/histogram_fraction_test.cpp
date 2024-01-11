// Copyright 2015-2018 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/fraction.hpp>
#include <boost/histogram/accumulators/ostream.hpp>
#include <boost/histogram/axis/integer.hpp>
#include <boost/histogram/histogram.hpp>
#include <boost/histogram/make_histogram.hpp>
#include <boost/histogram/ostream.hpp>
#include <valarray>
#include <vector>
#include "throw_exception.hpp"

using namespace boost::histogram;

int main() {
  using fraction_t = accumulators::fraction<>;

  auto h = make_histogram_with(dense_storage<fraction_t>(), axis::integer<>(0, 2));

  h(0, sample(true));
  h(0, sample(false));
  h(1, sample(true));

  BOOST_TEST_EQ(h.at(0), (fraction_t{1, 1}));
  BOOST_TEST_EQ(h.at(1), (fraction_t{1, 0}));

  // cannot use std::vector<bool> because of vector specialization
  std::vector<char> s = {true, false};
  std::vector<int> x = {0, 1};
  h.fill(x, sample(s));
  BOOST_TEST_EQ(h.at(0), (fraction_t{2, 1}));
  BOOST_TEST_EQ(h.at(1), (fraction_t{1, 1}));

  // any contiguous container of bool works which is not specialized
  std::array<bool, 2> s2 = {false, true};
  h.fill(x, sample(s2));
  BOOST_TEST_EQ(h.at(0), (fraction_t{2, 2}));
  BOOST_TEST_EQ(h.at(1), (fraction_t{2, 1}));

  std::valarray<bool> s3 = {false, true};
  h.fill(x, sample(s3));
  BOOST_TEST_EQ(h.at(0), (fraction_t{2, 3}));
  BOOST_TEST_EQ(h.at(1), (fraction_t{3, 1}));

  return boost::report_errors();
}
