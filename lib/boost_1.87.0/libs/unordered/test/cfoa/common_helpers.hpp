// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
  
#ifndef BOOST_UNORDERED_TEST_CFOA_COMMON_HELPERS_HPP
#define BOOST_UNORDERED_TEST_CFOA_COMMON_HELPERS_HPP

#include <boost/unordered/concurrent_flat_map_fwd.hpp>
#include <boost/unordered/concurrent_flat_set_fwd.hpp>
#include <boost/unordered/concurrent_node_map_fwd.hpp>
#include <boost/unordered/concurrent_node_set_fwd.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_node_map.hpp>
#include <boost/unordered/unordered_node_set.hpp>

#include <cstddef>
#include <type_traits>
#include <utility>

template <typename K>
struct value_cardinality
{
  static constexpr std::size_t value=1;
};

template <typename K, typename V>
struct value_cardinality<std::pair<K, V> >
{
  static constexpr std::size_t value=2;
};

template <typename K>
struct value_nonconst_cardinality
{
  static constexpr std::size_t value=1;
};

template <typename K, typename V>
struct value_nonconst_cardinality<std::pair<K, V> >
{
  static constexpr std::size_t value=
    1 * !std::is_const<K>::value + 
    1 * !std::is_const<V>::value ;
};

template <class Container>
struct is_container_node_based: std::false_type {};

template <typename K, typename V, typename H, typename P, typename A>
struct is_container_node_based<boost::concurrent_node_map<K, V, H, P, A> >
  : std::true_type {};

template <typename K, typename H, typename P, typename A>
struct is_container_node_based<boost::concurrent_node_set<K, H, P, A> >
  : std::true_type {};

template <class Container>
struct reference_container_impl;

template <class Container>
using reference_container = typename reference_container_impl<Container>::type;

template <typename K, typename V, typename H, typename P, typename A>
struct reference_container_impl<boost::concurrent_flat_map<K, V, H, P, A> >
{
  using type = boost::unordered_flat_map<K, V>;
};

template <typename K, typename V, typename H, typename P, typename A>
struct reference_container_impl<boost::concurrent_node_map<K, V, H, P, A> >
{
  using type = boost::unordered_node_map<K, V>;
};

template <typename K, typename H, typename P, typename A>
struct reference_container_impl<boost::concurrent_flat_set<K, H, P, A> >
{
  using type = boost::unordered_flat_set<K>;
};

template <typename K, typename H, typename P, typename A>
struct reference_container_impl<boost::concurrent_node_set<K, H, P, A> >
{
  using type = boost::unordered_node_set<K>;
};

template <class Container>
struct nonconcurrent_container_impl;

template <class Container>
using nonconcurrent_container = 
  typename nonconcurrent_container_impl<Container>::type;

template <typename K, typename V, typename H, typename P, typename A>
struct nonconcurrent_container_impl<boost::concurrent_flat_map<K, V, H, P, A> >
{
  using type = boost::unordered_flat_map<K, V, H, P, A>;
};

template <typename K, typename V, typename H, typename P, typename A>
struct nonconcurrent_container_impl<boost::concurrent_node_map<K, V, H, P, A> >
{
  using type = boost::unordered_node_map<K, V, H, P, A>;
};

template <typename K, typename H, typename P, typename A>
struct nonconcurrent_container_impl<boost::concurrent_flat_set<K, H, P, A> >
{
  using type = boost::unordered_flat_set<K, H, P, A>;
};

template <typename K, typename H, typename P, typename A>
struct nonconcurrent_container_impl<boost::concurrent_node_set<K, H, P, A> >
{
  using type = boost::unordered_node_set<K, H, P, A>;
};

template <typename Container, template <typename> class Allocator>
struct replace_allocator_impl;

template <typename Container, template <typename> class Allocator>
using replace_allocator = 
  typename replace_allocator_impl<Container, Allocator>::type;

template <
  typename K, typename V, typename H, typename P, typename A,
  template <typename> class Allocator
>
struct replace_allocator_impl<
  boost::concurrent_flat_map<K, V, H, P, A>, Allocator>
{
  using value_type = 
    typename boost::concurrent_flat_map<K, V, H, P, A>::value_type;
  using type = 
    boost::concurrent_flat_map<K, V, H, P, Allocator<value_type> >;
};

template <
  typename K, typename H, typename P, typename A,
  template <typename> class Allocator
>
struct replace_allocator_impl<
  boost::concurrent_flat_set<K, H, P, A>, Allocator>
{
  using value_type = 
    typename boost::concurrent_flat_set<K, H, P, A>::value_type;
  using type = 
    boost::concurrent_flat_set<K, H, P, Allocator<value_type> >;
};

template <
  typename K, typename V, typename H, typename P, typename A,
  template <typename> class Allocator
>
struct replace_allocator_impl<
  boost::concurrent_node_map<K, V, H, P, A>, Allocator>
{
  using value_type = 
    typename boost::concurrent_node_map<K, V, H, P, A>::value_type;
  using type = 
    boost::concurrent_node_map<K, V, H, P, Allocator<value_type> >;
};

template <
  typename K, typename H, typename P, typename A,
  template <typename> class Allocator
>
struct replace_allocator_impl<
  boost::concurrent_node_set<K, H, P, A>, Allocator>
{
  using value_type = 
    typename boost::concurrent_node_set<K, H, P, A>::value_type;
  using type = 
    boost::concurrent_node_set<K, H, P, Allocator<value_type> >;
};

template <typename K>
K const& get_key(K const& x) { return x; }

template <typename K,typename V>
K const& get_key(const std::pair<K, V>& x) { return x.first; }

template <typename K>
K const& get_value(K const& x) { return x; }

template <typename K,typename V>
V const& get_value(const std::pair<K, V>& x) { return x.second; }

template <typename K,typename V>
V& get_value(std::pair<K, V>& x) { return x.second; }

template <class X, class Y>
void test_matches_reference(X const& x, Y const& reference_cont)
{
  using value_type = typename X::value_type;
  BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& v) {
    BOOST_TEST(reference_cont.contains(get_key(v)));
    BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
  }));
}

template <class X, class Y>
void test_fuzzy_matches_reference(
  X const& x, Y const& reference_cont, test::random_generator rg)
{
  using value_type = typename X::value_type;
  BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& v) {
    BOOST_TEST(reference_cont.contains(get_key(v)));
    if (rg == test::sequential) {
      BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
    }
  }));
}

#endif // BOOST_UNORDERED_TEST_CFOA_COMMON_HELPERS_HPP
