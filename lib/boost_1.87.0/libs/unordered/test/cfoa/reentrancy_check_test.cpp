// Copyright 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <cstdlib>

#define BOOST_ENABLE_ASSERT_HANDLER 

static bool reentrancy_detected = false;

namespace boost {
  // Caveat lector: a proper handler shouldn't throw as it may be executed
  // within a noexcept function.

  void assertion_failed_msg(
    char const*, char const*, char const*, char const*, long)
  {
    reentrancy_detected = true;
    throw 0;
  }

  // LCOV_EXCL_START
  void assertion_failed(char const*, char const*, char const*, long)
  {
    std::abort();
  }
  // LCOV_EXCL_STOP

} // namespace boost

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>
#include <boost/core/lightweight_test.hpp>

using test::default_generator;

using map_type = boost::unordered::concurrent_flat_map<raii, raii>;
using node_map_type = boost::unordered::concurrent_node_map<raii, raii>;
using set_type = boost::unordered::concurrent_flat_set<raii>;
using node_set_type = boost::unordered::concurrent_node_set<raii>;

map_type* test_map;
map_type* test_node_map;
set_type* test_set;
node_set_type* test_node_set;

template<typename F>
void detect_reentrancy(F f)
{
  reentrancy_detected = false;
  try {
    f();
  }
  catch(int) {}
  BOOST_TEST(reentrancy_detected);
}

namespace {
  template <class X, class GF>
  void reentrancy_tests(X*, GF gen_factory, test::random_generator rg)
  {
    using key_type = typename X::key_type;

    // concurrent_flat_set visit is always const access
    using arg_type = typename std::conditional<
      std::is_same<typename X::key_type, typename X::value_type>::value,
      typename X::value_type const,
      typename X::value_type
    >::type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    X x1, x2;
    x1.insert(values.begin(), values.end());
    x2.insert(values.begin(), values.end());

    detect_reentrancy([&] {
      x1.visit_all([&](arg_type&) { (void)x1.contains(key_type()); });
    }); // LCOV_EXCL_LINE

    detect_reentrancy([&] {
      x1.visit_all([&](arg_type&) { x1.rehash(0); });
    }); // LCOV_EXCL_LINE

    detect_reentrancy([&] {
      x1.visit_all([&](arg_type&) { 
        x2.visit_all([&](arg_type&) { 
          x1=x2;
        }); // LCOV_EXCL_START
      });
    });
    // LCOV_EXCL_STOP

    detect_reentrancy([&] {
      x1.visit_all([&](arg_type&) { 
        x2.visit_all([&](arg_type&) { 
          x2=x1;
        }); // LCOV_EXCL_START
      });
    });
    // LCOV_EXCL_STOP
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  reentrancy_tests,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)))
// clang-format on

RUN_TESTS()
