
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
template <typename T, typename H, typename P, typename A>
class instantiate_flat_set
{
  typedef boost::unordered_flat_set<T, H, P, A> container;
  container x;
};

template class instantiate_flat_set<int, boost::hash<int>, std::equal_to<int>,
  test::minimal::allocator<int> >;

template class instantiate_flat_set<test::minimal::assignable const,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;

template <typename T, typename H, typename P, typename A>
class instantiate_node_set
{
  typedef boost::unordered_node_set<T, H, P, A> container;
  container x;
};

template class instantiate_node_set<int, boost::hash<int>, std::equal_to<int>,
  test::minimal::allocator<int> >;

template class instantiate_node_set<test::minimal::assignable const,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;

#else

#define INSTANTIATE(type)                                                      \
  template class boost::unordered::detail::instantiate_##type

INSTANTIATE(set)<int, boost::hash<int>, std::equal_to<int>,
  test::minimal::allocator<int> >;
INSTANTIATE(multiset)<int const, boost::hash<int>, std::equal_to<int>,
  test::minimal::allocator<int> >;

INSTANTIATE(set)<test::minimal::assignable const,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;
INSTANTIATE(multiset)<test::minimal::assignable,
  test::minimal::hash<test::minimal::assignable>,
  test::minimal::equal_to<test::minimal::assignable>,
  test::minimal::allocator<int> >;

#endif

template <class X> static void type_traits_impl()
{
  BOOST_STATIC_ASSERT(boost::is_same<int const&,
    typename std::iterator_traits<typename X::iterator>::reference>::value);
}

UNORDERED_AUTO_TEST (type_traits) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  type_traits_impl<boost::unordered_flat_set<int> >();
  type_traits_impl<boost::unordered_node_set<int> >();
#else
  type_traits_impl<boost::unordered_set<int> >();
  type_traits_impl<boost::unordered_multiset<int> >();
#endif
}

template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void test0_impl()
{
  test::minimal::constructor_param x;

  test::minimal::assignable assignable(x);

  Set<int> int_set;

  Set<int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<int> >
    int_set2;

  Set<test::minimal::assignable, test::minimal::hash<test::minimal::assignable>,
    test::minimal::equal_to<test::minimal::assignable>,
    test::minimal::allocator<test::minimal::assignable> >
    set;

  container_test(int_set, 0);
  container_test(int_set2, 0);
  container_test(set, assignable);
}

UNORDERED_AUTO_TEST (test0) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  test0_impl<boost::unordered_flat_set>();
  test0_impl<boost::unordered_node_set>();
#else
  test0_impl<boost::unordered_set>();
  test0_impl<boost::unordered_multiset>();
#endif
}

template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void equality_tests_impl()
{
  typedef test::minimal::copy_constructible_equality_comparable value_type;

  Set<int> int_set;

  Set<int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<int> >
    int_set2;

  Set<test::minimal::copy_constructible_equality_comparable,
    test::minimal::hash<test::minimal::copy_constructible_equality_comparable>,
    test::minimal::equal_to<
      test::minimal::copy_constructible_equality_comparable>,
    test::minimal::allocator<value_type> >
    set;

  equality_test(int_set);
  equality_test(int_set2);
  equality_test(set);
}

UNORDERED_AUTO_TEST (equality_tests) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  equality_tests_impl<boost::unordered_flat_set>();
  equality_tests_impl<boost::unordered_node_set>();
#else
  equality_tests_impl<boost::unordered_set>();
  equality_tests_impl<boost::unordered_multiset>();
#endif
}

template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void test1_unique_impl()
{

  boost::hash<int> hash;
  std::equal_to<int> equal_to;
  int value = 0;
  Set<int> set;

  Set<int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<int> >
    set2;

  unordered_unique_test(set, value);
  unordered_set_test(set, value);
  unordered_copyable_test(set, value, value, hash, equal_to);

  unordered_unique_test(set2, value);
  unordered_set_test(set2, value);
  unordered_copyable_test(set2, value, value, hash, equal_to);
}

#ifndef BOOST_UNORDERED_FOA_TESTS
template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void test1_equivalent_impl()
{

  boost::hash<int> hash;
  std::equal_to<int> equal_to;
  int value = 0;
  Set<int> set;

  Set<int, boost::hash<int>, std::equal_to<int>,
    test::minimal::cxx11_allocator<int> >
    set2;

  unordered_equivalent_test(set, value);
  unordered_set_test(set, value);
  unordered_copyable_test(set, value, value, hash, equal_to);

  unordered_equivalent_test(set2, value);
  unordered_set_test(set2, value);
  unordered_copyable_test(set2, value, value, hash, equal_to);
}
#endif

UNORDERED_AUTO_TEST (test1) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  test1_unique_impl<boost::unordered_flat_set>();
  test1_unique_impl<boost::unordered_node_set>();
#else
  test1_unique_impl<boost::unordered_set>();
  test1_equivalent_impl<boost::unordered_multiset>();
#endif
}

template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void test2_unique_impl()
{
  test::minimal::constructor_param x;

  test::minimal::assignable assignable(x);
  test::minimal::copy_constructible copy_constructible(x);
  test::minimal::hash<test::minimal::assignable> hash(x);
  test::minimal::equal_to<test::minimal::assignable> equal_to(x);

  Set<test::minimal::assignable, test::minimal::hash<test::minimal::assignable>,
    test::minimal::equal_to<test::minimal::assignable>,
    test::minimal::allocator<test::minimal::assignable> >
    set;

  unordered_unique_test(set, assignable);
  unordered_set_test(set, assignable);
  unordered_copyable_test(set, assignable, assignable, hash, equal_to);
  unordered_set_member_test(set, assignable);
}

#ifndef BOOST_UNORDERED_FOA_TESTS
template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void test2_equivalent_impl()
{
  test::minimal::constructor_param x;

  test::minimal::assignable assignable(x);
  test::minimal::copy_constructible copy_constructible(x);
  test::minimal::hash<test::minimal::assignable> hash(x);
  test::minimal::equal_to<test::minimal::assignable> equal_to(x);

  Set<test::minimal::assignable, test::minimal::hash<test::minimal::assignable>,
    test::minimal::equal_to<test::minimal::assignable>,
    test::minimal::allocator<test::minimal::assignable> >
    set;

  unordered_equivalent_test(set, assignable);
  unordered_set_test(set, assignable);
  unordered_copyable_test(set, assignable, assignable, hash, equal_to);
  unordered_set_member_test(set, assignable);
}
#endif

UNORDERED_AUTO_TEST (test2) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  test2_unique_impl<boost::unordered_flat_set>();
  test2_unique_impl<boost::unordered_node_set>();
#else
  test2_unique_impl<boost::unordered_set>();
  test2_equivalent_impl<boost::unordered_multiset>();
#endif
}

template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void movable1_tests_impl()
{
  test::minimal::constructor_param x;

  test::minimal::movable1 movable1(x);
  test::minimal::hash<test::minimal::movable1> hash(x);
  test::minimal::equal_to<test::minimal::movable1> equal_to(x);

  Set<test::minimal::movable1, test::minimal::hash<test::minimal::movable1>,
    test::minimal::equal_to<test::minimal::movable1>,
    test::minimal::allocator<test::minimal::movable1> >
    set;

  // TODO: find out why Daniel had this commented out and if we need it and the
  // corresponding equivalent impl
  //
  // unordered_unique_test(set, movable1);
  unordered_set_test(set, movable1);
  unordered_movable_test(set, movable1, movable1, hash, equal_to);
}

UNORDERED_AUTO_TEST (movable1_tests) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  movable1_tests_impl<boost::unordered_flat_set>();
  movable1_tests_impl<boost::unordered_node_set>();
#else
  movable1_tests_impl<boost::unordered_set>();
  movable1_tests_impl<boost::unordered_multiset>();
#endif
}

template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void movable2_tests_impl()
{
  test::minimal::constructor_param x;

  test::minimal::movable2 movable2(x);
  test::minimal::hash<test::minimal::movable2> hash(x);
  test::minimal::equal_to<test::minimal::movable2> equal_to(x);

  Set<test::minimal::movable2, test::minimal::hash<test::minimal::movable2>,
    test::minimal::equal_to<test::minimal::movable2>,
    test::minimal::allocator<test::minimal::movable2> >
    set;

  // unordered_unique_test(set, movable2);
  unordered_set_test(set, movable2);
  unordered_movable_test(set, movable2, movable2, hash, equal_to);
}

UNORDERED_AUTO_TEST (movable2_tests) {

#ifdef BOOST_UNORDERED_FOA_TESTS
  movable2_tests_impl<boost::unordered_flat_set>();
  movable2_tests_impl<boost::unordered_node_set>();
#else
  movable2_tests_impl<boost::unordered_set>();
  movable2_tests_impl<boost::unordered_multiset>();
#endif
}

template <template <class T, class H = boost::hash<T>,
  class P = std::equal_to<T>, class A = std::allocator<T> >
  class Set>
static void destructible_tests_impl()
{
  test::minimal::constructor_param x;

  test::minimal::destructible destructible(x);
  test::minimal::hash<test::minimal::destructible> hash(x);
  test::minimal::equal_to<test::minimal::destructible> equal_to(x);

  Set<test::minimal::destructible,
    test::minimal::hash<test::minimal::destructible>,
    test::minimal::equal_to<test::minimal::destructible> >
    set;

  unordered_destructible_test(set);
}

UNORDERED_AUTO_TEST (destructible_tests) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  destructible_tests_impl<boost::unordered_flat_set>();
  destructible_tests_impl<boost::unordered_node_set>();
#else
  destructible_tests_impl<boost::unordered_set>();
  destructible_tests_impl<boost::unordered_multiset>();
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
    boost::unordered_flat_set<lwg2059_key> x;
    x.emplace(lwg2059_key(10));
    x.erase(x.begin());
  }

  {
    boost::unordered_node_set<lwg2059_key> x;
    x.emplace(lwg2059_key(10));
    x.erase(x.begin());
  }
#else
  {
    boost::unordered_set<lwg2059_key> x;
    x.emplace(lwg2059_key(10));
    x.erase(x.begin());
  }

  {
    boost::unordered_multiset<lwg2059_key> x;
    x.emplace(lwg2059_key(10));
    x.erase(x.begin());
  }
#endif
}

RUN_TESTS()
