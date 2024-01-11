// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"
#include <boost/config/workaround.hpp>
#include <boost/unordered/concurrent_flat_map_fwd.hpp>
#include <boost/unordered/concurrent_flat_set_fwd.hpp>
#include <limits>

test::seed_t initialize_seed{32304628};

using test::default_generator;
using test::limited_range;
using test::sequential;

template <class T>
void swap_call(boost::unordered::concurrent_flat_map<T, T>& x1,
  boost::unordered::concurrent_flat_map<T, T>& x2)
{
  swap(x1, x2);
}

template <class T>
bool equal_call(boost::unordered::concurrent_flat_map<T, T>& x1,
  boost::unordered::concurrent_flat_map<T, T>& x2)
{
  return x1 == x2;
}

template <class T>
bool unequal_call(boost::unordered::concurrent_flat_map<T, T>& x1,
  boost::unordered::concurrent_flat_map<T, T>& x2)
{
  return x1 != x2;
}

template <class T>
void swap_call(boost::unordered::concurrent_flat_set<T>& x1,
  boost::unordered::concurrent_flat_set<T>& x2)
{
  swap(x1, x2);
}

template <class T>
bool equal_call(boost::unordered::concurrent_flat_set<T>& x1,
  boost::unordered::concurrent_flat_set<T>& x2)
{
  return x1 == x2;
}

template <class T>
bool unequal_call(boost::unordered::concurrent_flat_set<T>& x1,
  boost::unordered::concurrent_flat_set<T>& x2)
{
  return x1 != x2;
}

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>

using map_type = boost::unordered::concurrent_flat_map<int, int>;
using set_type = boost::unordered::concurrent_flat_map<int, int>;

map_type* test_map;
set_type* test_set;

template <typename X>
void fwd_swap_call(X*)
{
#if !defined(BOOST_CLANG_VERSION) ||                                          \
  BOOST_WORKAROUND(BOOST_CLANG_VERSION, < 30700) ||                           \
  BOOST_WORKAROUND(BOOST_CLANG_VERSION, >= 30800)
// clang-3.7 seems to have a codegen bug here so we workaround it

  X x1, x2;
  swap_call(x1, x2);
#endif
}

template <typename X>
void fwd_equal_call(X*)
{
  X x1, x2;
  BOOST_TEST(equal_call(x1, x2));
}

template <typename X>
void fwd_unequal_call(X*)
{
  X x1, x2;
  BOOST_TEST_NOT(unequal_call(x1, x2));
}

// this isn't the best place for this test but it's better than introducing a
// new file
template <typename X>
void max_size(X*)
{
  X x1;
  BOOST_TEST_EQ(
    x1.max_size(), std::numeric_limits<typename X::size_type>::max());
}

// clang-format off
UNORDERED_TEST(
  fwd_swap_call,
  ((test_map)(test_set)))

UNORDERED_TEST(
  fwd_equal_call,
  ((test_map)(test_set)))

UNORDERED_TEST(
  fwd_unequal_call,
  ((test_map)(test_set)))

UNORDERED_TEST(
  max_size,
  ((test_map)(test_set)))
// clang-format on

RUN_TESTS()
