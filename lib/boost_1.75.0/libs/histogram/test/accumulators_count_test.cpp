// Copyright 2015-2018 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/count.hpp>
#include <boost/histogram/accumulators/ostream.hpp>
#include "throw_exception.hpp"
#include "utility_str.hpp"

using namespace boost::histogram;
using namespace std::literals;

int main() {
  using c_t = accumulators::count<double>;

  c_t c;
  ++c;
  BOOST_TEST_EQ(c.value(), 1);
  BOOST_TEST_EQ(str(c), "1"s);
  BOOST_TEST_EQ(str(c, 2, false), " 1"s);
  BOOST_TEST_EQ(str(c, 2, true), "1 "s);

  c += 2;
  BOOST_TEST_EQ(str(c), "3"s);

  BOOST_TEST_EQ(c, 3);
  BOOST_TEST_NE(c, 2);

  c_t one(1), two(2), one_copy(1);
  BOOST_TEST_LT(one, two);
  BOOST_TEST_LE(one, two);
  BOOST_TEST_LE(one, one_copy);
  BOOST_TEST_GT(two, one);
  BOOST_TEST_GE(two, one);
  BOOST_TEST_GE(one, one_copy);

  BOOST_TEST_EQ(c_t{} += c_t{}, c_t{});

  return boost::report_errors();
}
