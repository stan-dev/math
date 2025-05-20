
// Copyright 2008-2009 Daniel James.
// Copyright 2022-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// clang-format off
#include "../helpers/prefix.hpp"
#ifdef BOOST_UNORDERED_FOA_TESTS
#include <boost/unordered/unordered_flat_set_fwd.hpp>
#include <boost/unordered/unordered_node_set_fwd.hpp>
#include <boost/unordered/detail/implementation.hpp>
#else
#include <boost/unordered/unordered_set_fwd.hpp>
#endif
#include "../helpers/postfix.hpp"
// clang-format on

struct true_type
{
  char x[100];
};
struct false_type
{
  char x;
};

false_type is_unordered_set_impl(void*);

#ifdef BOOST_UNORDERED_FOA_TESTS
template <class Value, class Hash, class Pred, class Alloc>
true_type is_unordered_set_impl(
  boost::unordered_flat_set<Value, Hash, Pred, Alloc>*);

template <class Value, class Hash, class Pred, class Alloc>
true_type is_unordered_set_impl(
  boost::unordered_node_set<Value, Hash, Pred, Alloc>*);

template <typename T>
void call_swap(boost::unordered_flat_set<T>& x, boost::unordered_flat_set<T>& y)
{
  swap(x, y);
}

template <typename T>
void call_swap(boost::unordered_node_set<T>& x, boost::unordered_node_set<T>& y)
{
  swap(x, y);
}

template <typename T>
bool call_equals(
  boost::unordered_flat_set<T>& x, boost::unordered_flat_set<T>& y)
{
  return x == y;
}

template <typename T>
bool call_equals(
  boost::unordered_node_set<T>& x, boost::unordered_node_set<T>& y)
{
  return x == y;
}

template <typename T>
bool call_not_equals(
  boost::unordered_flat_set<T>& x, boost::unordered_flat_set<T>& y)
{
  return x != y;
}

template <typename T>
bool call_not_equals(
  boost::unordered_node_set<T>& x, boost::unordered_node_set<T>& y)
{
  return x != y;
}
#else
template <class Value, class Hash, class Pred, class Alloc>
true_type is_unordered_set_impl(
  boost::unordered_set<Value, Hash, Pred, Alloc>*);

template <typename T>
void call_swap(boost::unordered_set<T>& x, boost::unordered_set<T>& y)
{
  swap(x, y);
}

template <typename T>
bool call_equals(boost::unordered_set<T>& x, boost::unordered_set<T>& y)
{
  return x == y;
}

template <typename T>
bool call_not_equals(boost::unordered_set<T>& x, boost::unordered_set<T>& y)
{
  return x != y;
}
#endif

#ifndef BOOST_UNORDERED_FOA_TESTS
template <typename T>
void call_swap(boost::unordered_multiset<T>& x, boost::unordered_multiset<T>& y)
{
  swap(x, y);
}

template <typename T>
bool call_equals(
  boost::unordered_multiset<T>& x, boost::unordered_multiset<T>& y)
{
  return x == y;
}

template <typename T>
bool call_not_equals(
  boost::unordered_multiset<T>& x, boost::unordered_multiset<T>& y)
{
  return x != y;
}
#endif

#include "../helpers/test.hpp"

#ifdef BOOST_UNORDERED_FOA_TESTS
typedef boost::unordered_flat_set<int> int_set;
typedef boost::unordered_node_set<int> int_node_set;
#else
typedef boost::unordered_set<int> int_set;
typedef boost::unordered_multiset<int> int_multiset;
#endif

UNORDERED_AUTO_TEST (use_fwd_declared_trait_without_definition) {
  BOOST_TEST(sizeof(is_unordered_set_impl((int_set*)0)) == sizeof(true_type));
#ifdef BOOST_UNORDERED_FOA_TESTS
  BOOST_TEST(
    sizeof(is_unordered_set_impl((int_node_set*)0)) == sizeof(true_type));
#endif
}

#ifdef BOOST_UNORDERED_FOA_TESTS
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_node_set.hpp>
#else
#include <boost/unordered_set.hpp>
#endif

UNORDERED_AUTO_TEST (use_fwd_declared_trait) {
  int_set x;
  BOOST_TEST(sizeof(is_unordered_set_impl(&x)) == sizeof(true_type));
  BOOST_TEST(sizeof(is_unordered_set_impl((int*)0)) == sizeof(false_type));
}

#ifdef BOOST_UNORDERED_FOA_TESTS
UNORDERED_AUTO_TEST (use_node_fwd_declared_trait) {
  int_node_set x;
  BOOST_TEST(sizeof(is_unordered_set_impl(&x)) == sizeof(true_type));
  BOOST_TEST(sizeof(is_unordered_set_impl((int*)0)) == sizeof(false_type));
}
#endif

UNORDERED_AUTO_TEST (use_set_fwd_declared_function) {
  int_set x, y;
  x.insert(1);
  y.insert(2);
  call_swap(x, y);

  BOOST_TEST(y.find(1) != y.end());
  BOOST_TEST(y.find(2) == y.end());

  BOOST_TEST(x.find(1) == x.end());
  BOOST_TEST(x.find(2) != x.end());

  BOOST_TEST(!call_equals(x, y));
  BOOST_TEST(call_not_equals(x, y));
}

#ifdef BOOST_UNORDERED_FOA_TESTS
UNORDERED_AUTO_TEST (use_node_set_fwd_declared_function) {
  int_node_set x, y;
  x.insert(1);
  y.insert(2);
  call_swap(x, y);

  BOOST_TEST(y.find(1) != y.end());
  BOOST_TEST(y.find(2) == y.end());

  BOOST_TEST(x.find(1) == x.end());
  BOOST_TEST(x.find(2) != x.end());

  BOOST_TEST(!call_equals(x, y));
  BOOST_TEST(call_not_equals(x, y));
}
#endif

#ifndef BOOST_UNORDERED_FOA_TESTS
UNORDERED_AUTO_TEST (use_multiset_fwd_declared_function) {
  int_multiset x, y;
  call_swap(x, y);
  BOOST_TEST(call_equals(x, y));
  BOOST_TEST(!call_not_equals(x, y));
}
#endif

RUN_TESTS()
