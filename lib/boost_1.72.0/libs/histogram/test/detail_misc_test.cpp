// Copyright 2015-2019 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/histogram/accumulators/weighted_sum.hpp>
#include <boost/histogram/detail/common_type.hpp>
#include <boost/histogram/detail/counting_streambuf.hpp>
#include <boost/histogram/detail/non_member_container_access.hpp>
#include <boost/histogram/fwd.hpp>
#include <boost/histogram/literals.hpp>
#include <boost/histogram/storage_adaptor.hpp>
#include <boost/histogram/unlimited_storage.hpp>
#include <ostream>
#include "std_ostream.hpp"

using namespace boost::histogram;
using namespace boost::histogram::literals;
namespace dtl = boost::histogram::detail;

int main() {
  // literals
  {
    BOOST_TEST_TRAIT_SAME(std::integral_constant<unsigned, 0>, decltype(0_c));
    BOOST_TEST_TRAIT_SAME(std::integral_constant<unsigned, 3>, decltype(3_c));
    BOOST_TEST_EQ(decltype(10_c)::value, 10);
    BOOST_TEST_EQ(decltype(213_c)::value, 213);
  }

  // common_storage
  {
    BOOST_TEST_TRAIT_SAME(dtl::common_storage<unlimited_storage<>, unlimited_storage<>>,
                          unlimited_storage<>);
    BOOST_TEST_TRAIT_SAME(
        dtl::common_storage<dense_storage<double>, dense_storage<double>>,
        dense_storage<double>);
    BOOST_TEST_TRAIT_SAME(dtl::common_storage<dense_storage<int>, dense_storage<double>>,
                          dense_storage<double>);
    BOOST_TEST_TRAIT_SAME(dtl::common_storage<dense_storage<double>, dense_storage<int>>,
                          dense_storage<double>);
    BOOST_TEST_TRAIT_SAME(dtl::common_storage<dense_storage<double>, unlimited_storage<>>,
                          dense_storage<double>);
    BOOST_TEST_TRAIT_SAME(dtl::common_storage<dense_storage<int>, unlimited_storage<>>,
                          unlimited_storage<>);
    BOOST_TEST_TRAIT_SAME(dtl::common_storage<dense_storage<double>, weight_storage>,
                          weight_storage);
  }

  // size & data
  {
    char a[4] = {1, 2, 3, 4};
    BOOST_TEST_EQ(dtl::size(a), 4u);
    BOOST_TEST_EQ(dtl::data(a), a);
    auto b = {1, 2};
    BOOST_TEST_EQ(dtl::size(b), 2u);
    BOOST_TEST_EQ(dtl::data(b), b.begin());
    struct C {
      unsigned size() const { return 3; }
      int* data() { return buf; }
      const int* data() const { return buf; }
      int buf[1];
    } c;
    BOOST_TEST_EQ(dtl::size(c), 3u);
    BOOST_TEST_EQ(dtl::data(c), c.buf);
    BOOST_TEST_EQ(dtl::data(static_cast<const C&>(c)), c.buf);
    struct {
      int size() const { return 5; }
    } d;
    BOOST_TEST_EQ(dtl::size(d), 5u);
  }

  // counting_streambuf
  {
    dtl::counting_streambuf<char> cbuf;
    std::ostream os(&cbuf);
    os.put('x');
    BOOST_TEST_EQ(cbuf.count, 1);
    os << 12;
    BOOST_TEST_EQ(cbuf.count, 3);
    os << "123";
    BOOST_TEST_EQ(cbuf.count, 6);
  }

  return boost::report_errors();
}
