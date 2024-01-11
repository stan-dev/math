
// Copyright 2006-2009 Daniel James.
// Copyright 2022-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// This test creates the containers with members that meet their minimum
// requirements. Makes sure everything compiles and is defined correctly.

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"
#include "../objects/minimal.hpp"
#include "./compile_tests.hpp"

// Explicit instantiation to catch compile-time errors
#ifdef BOOST_UNORDERED_FOA_TESTS

// emulates what was already done for previous tests but without leaking to
// the detail namespace
//
template <typename K, typename T, typename H, typename P, typename A>
class instantiate_flat_map
{
  typedef boost::unordered_flat_map<K, T, H, P, A> container;
  container x;
};

template class instantiate_flat_map<int, int, boost::hash<int>,
  std::equal_to<int>, test::minimal::allocator<int> >;

template class instantiate_flat_map<test::minimal::assignable const,
  test::minimal::default_assignable const,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;

template <typename K, typename T, typename H, typename P, typename A>
class instantiate_node_map
{
  typedef boost::unordered_node_map<K, T, H, P, A> container;
  container x;
};

template class instantiate_node_map<int, int, boost::hash<int>,
  std::equal_to<int>, test::minimal::allocator<int> >;

template class instantiate_node_map<test::minimal::assignable const,
  test::minimal::default_assignable const,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;

#else
#define INSTANTIATE(type)                                                      \
  template class boost::unordered::detail::instantiate_##type

INSTANTIATE(map)<int, int, boost::hash<int>, std::equal_to<int>,
  test::minimal::allocator<int> >;
INSTANTIATE(multimap)<int const, int const, boost::hash<int>,
  std::equal_to<int>, test::minimal::allocator<int> >;

INSTANTIATE(
  map)<test::minimal::assignable const, test::minimal::default_assignable const,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;
INSTANTIATE(multimap)<test::minimal::assignable, test::minimal::assignable,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;
#endif

template <template <class Key, class T, class H = boost::hash<Key>,
  class P = std::equal_to<Key>,
  class Allocator = std::allocator<std::pair<Key const, T> > >
  class Map>
static void test0_impl()
{
  test::minimal::constructor_param x;

  typedef std::pair<test::minimal::assignable const, test::minimal::assignable>
    value_type;
  value_type value(x, x);

  Map<int, int> int_map;

  Map<int, int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<std::pair<int const, int> > >
    int_map2;

  Map<test::minimal::assignable, test::minimal::assignable,
    test::minimal::hash<test::minimal::assignable>,
    test::minimal::equal_to<test::minimal::assignable>,
    test::minimal::allocator<value_type> >
    map;

  container_test(int_map, std::pair<int const, int>(0, 0));
  container_test(int_map2, std::pair<int const, int>(0, 0));
  container_test(map, value);
}

UNORDERED_AUTO_TEST (test0) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  test0_impl<boost::unordered_flat_map>();
  test0_impl<boost::unordered_node_map>();
#else
  test0_impl<boost::unordered_map>();
  test0_impl<boost::unordered_multimap>();
#endif
}

template <template <class Key, class T, class H = boost::hash<Key>,
  class P = std::equal_to<Key>,
  class Allocator = std::allocator<std::pair<Key const, T> > >
  class Map>
static void equality_tests_impl()
{
  typedef std::pair<test::minimal::copy_constructible_equality_comparable const,
    test::minimal::copy_constructible_equality_comparable>
    value_type;

  Map<int, int> int_map;

  Map<int, int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<std::pair<int const, int> > >
    int_map2;

  Map<test::minimal::copy_constructible_equality_comparable,
    test::minimal::copy_constructible_equality_comparable,
    test::minimal::hash<test::minimal::copy_constructible_equality_comparable>,
    test::minimal::equal_to<
      test::minimal::copy_constructible_equality_comparable>,
    test::minimal::allocator<value_type> >
    map;

  equality_test(int_map);
  equality_test(int_map2);
  equality_test(map);
}

UNORDERED_AUTO_TEST (equality_tests) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  equality_tests_impl<boost::unordered_flat_map>();
  equality_tests_impl<boost::unordered_node_map>();
#else
  equality_tests_impl<boost::unordered_map>();
  equality_tests_impl<boost::unordered_multimap>();
#endif
}

template <template <class Key, class T, class H = boost::hash<Key>,
  class P = std::equal_to<Key>,
  class Allocator = std::allocator<std::pair<Key const, T> > >
  class Map>
static void test1_unique_impl()
{

  boost::hash<int> hash;
  std::equal_to<int> equal_to;
  int value = 0;
  std::pair<int const, int> map_value(0, 0);

  Map<int, int> map;

  Map<int, int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<std::pair<int const, int> > >
    map2;

  unordered_unique_test(map, map_value);
  unordered_map_test(map, value, value);
  unordered_copyable_test(map, value, map_value, hash, equal_to);
  unordered_map_functions(map, value, value);

  unordered_unique_test(map2, map_value);
  unordered_map_test(map2, value, value);
  unordered_copyable_test(map2, value, map_value, hash, equal_to);
  unordered_map_functions(map2, value, value);
}

#ifndef BOOST_UNORDERED_FOA_TESTS
template <template <class Key, class T, class H = boost::hash<Key>,
  class P = std::equal_to<Key>,
  class Allocator = std::allocator<std::pair<Key const, T> > >
  class Map>
static void test1_equivalent_impl()
{

  boost::hash<int> hash;
  std::equal_to<int> equal_to;
  int value = 0;
  std::pair<int const, int> map_value(0, 0);

  Map<int, int> map;

  Map<int, int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<std::pair<int const, int> > >
    map2;

  unordered_equivalent_test(map, map_value);
  unordered_map_test(map, value, value);
  unordered_copyable_test(map, value, map_value, hash, equal_to);

  unordered_equivalent_test(map2, map_value);
  unordered_map_test(map2, value, value);
  unordered_copyable_test(map2, value, map_value, hash, equal_to);
}
#endif

UNORDERED_AUTO_TEST (test1) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  test1_unique_impl<boost::unordered_flat_map>();
  test1_unique_impl<boost::unordered_node_map>();
#else
  test1_unique_impl<boost::unordered_map>();
  test1_equivalent_impl<boost::unordered_multimap>();
#endif
}

template <template <class Key, class T, class H = boost::hash<Key>,
  class P = std::equal_to<Key>,
  class Allocator = std::allocator<std::pair<Key const, T> > >
  class Map>
static void test2_unique_impl()
{
  test::minimal::constructor_param x;

  test::minimal::assignable assignable(x);
  test::minimal::copy_constructible copy_constructible(x);
  test::minimal::hash<test::minimal::assignable> hash(x);
  test::minimal::equal_to<test::minimal::assignable> equal_to(x);

  typedef std::pair<test::minimal::assignable const, test::minimal::assignable>
    map_value_type;
  map_value_type map_value(assignable, assignable);

  Map<test::minimal::assignable, test::minimal::assignable,
    test::minimal::hash<test::minimal::assignable>,
    test::minimal::equal_to<test::minimal::assignable>,
    test::minimal::allocator<map_value_type> >
    map;

  unordered_unique_test(map, map_value);
  unordered_map_test(map, assignable, assignable);
  unordered_copyable_test(map, assignable, map_value, hash, equal_to);
  unordered_map_member_test(map, map_value);

  Map<test::minimal::assignable, test::minimal::default_assignable,
    test::minimal::hash<test::minimal::assignable>,
    test::minimal::equal_to<test::minimal::assignable>,
    test::minimal::allocator<map_value_type> >
    map2;

  test::minimal::default_assignable default_assignable;

  unordered_map_functions(map2, assignable, default_assignable);
}

#ifndef BOOST_UNORDERED_FOA_TESTS
template <template <class Key, class T, class H = boost::hash<Key>,
  class P = std::equal_to<Key>,
  class Allocator = std::allocator<std::pair<Key const, T> > >
  class Map>
static void test2_equivalent_impl()
{
  test::minimal::constructor_param x;

  test::minimal::assignable assignable(x);
  test::minimal::copy_constructible copy_constructible(x);
  test::minimal::hash<test::minimal::assignable> hash(x);
  test::minimal::equal_to<test::minimal::assignable> equal_to(x);

  typedef std::pair<test::minimal::assignable const, test::minimal::assignable>
    map_value_type;
  map_value_type map_value(assignable, assignable);

  Map<test::minimal::assignable, test::minimal::assignable,
    test::minimal::hash<test::minimal::assignable>,
    test::minimal::equal_to<test::minimal::assignable>,
    test::minimal::allocator<map_value_type> >
    map;

  unordered_equivalent_test(map, map_value);
  unordered_map_test(map, assignable, assignable);
  unordered_copyable_test(map, assignable, map_value, hash, equal_to);
  unordered_map_member_test(map, map_value);
}
#endif

UNORDERED_AUTO_TEST (test2) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  test2_unique_impl<boost::unordered_flat_map>();
  test2_unique_impl<boost::unordered_node_map>();
#else
  test2_unique_impl<boost::unordered_map>();
  test2_equivalent_impl<boost::unordered_multimap>();
#endif
}

// Test for ambiguity when using key convertible from iterator
// See LWG2059

struct lwg2059_key
{
  int value;

  template <typename T> lwg2059_key(T v) : value(v) {}
};

std::size_t hash_value(lwg2059_key x)
{
  return static_cast<std::size_t>(x.value);
}

bool operator==(lwg2059_key x, lwg2059_key y) { return x.value == y.value; }

UNORDERED_AUTO_TEST (lwg2059) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  {
    boost::unordered_flat_map<lwg2059_key, int> x;
    x.emplace(lwg2059_key(10), 5);
    x.erase(x.begin());
  }

  {
    boost::unordered_node_map<lwg2059_key, int> x;
    x.emplace(lwg2059_key(10), 5);
    x.erase(x.begin());
  }
#else
  {
    boost::unordered_map<lwg2059_key, int> x;
    x.emplace(lwg2059_key(10), 5);
    x.erase(x.begin());
  }

  {
    boost::unordered_multimap<lwg2059_key, int> x;
    x.emplace(lwg2059_key(10), 5);
    x.erase(x.begin());
  }
#endif
}

RUN_TESTS()
