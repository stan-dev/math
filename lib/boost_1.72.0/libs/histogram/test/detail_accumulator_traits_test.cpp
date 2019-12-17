// Copyright 2019 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/histogram/detail/accumulator_traits.hpp>
#include <boost/histogram/weight.hpp>
#include <tuple>

namespace dtl = boost::histogram::detail;

int main() {
  using boost::histogram::weight_type;

  struct A1 {
    void operator()(){};
  };

  BOOST_TEST_NOT(dtl::accumulator_traits<A1>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<A1>::args, std::tuple<>);

  struct A2 {
    void operator()(int, double) {}
  };

  BOOST_TEST_NOT(dtl::accumulator_traits<A2>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<A2>::args,
                        std::tuple<int, double>);

  struct A3 {
    void operator()() {}
    void operator()(weight_type<int>) {}
  };

  BOOST_TEST(dtl::accumulator_traits<A3>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<A3>::args, std::tuple<>);

  struct A4 {
    void operator()(int, double, char) {}
    void operator()(weight_type<int>, int, double, char) {}
  };

  BOOST_TEST(dtl::accumulator_traits<A4>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<A4>::args,
                        std::tuple<int, double, char>);

  struct A5 {
    void operator()(const weight_type<int>, int) {}
  };

  BOOST_TEST(dtl::accumulator_traits<A5>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<A5>::args, std::tuple<int>);

  struct A6 {
    void operator()(const weight_type<int>&, const int) {}
  };

  BOOST_TEST(dtl::accumulator_traits<A6>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<A6>::args, std::tuple<int>);

  struct A7 {
    void operator()(weight_type<int>&&, int&&) {}
  };

  BOOST_TEST(dtl::accumulator_traits<A7>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<A7>::args, std::tuple<int&&>);

  struct B {
    int operator+=(int) { return 0; }
  };

  BOOST_TEST(dtl::accumulator_traits<B>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(typename dtl::accumulator_traits<B>::args, std::tuple<>);

  BOOST_TEST(dtl::accumulator_traits<int>::wsupport::value);
  BOOST_TEST_TRAIT_SAME(dtl::accumulator_traits<int>::args, std::tuple<>);

  return boost::report_errors();
}
