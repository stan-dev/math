
// Copyright 2024 Joaquin M Lopez Munz.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or move at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNORDERED_TEST_REPLACE_ALLOCATOR
#define BOOST_UNORDERED_TEST_REPLACE_ALLOCATOR

#include <boost/core/allocator_access.hpp> 
#include <utility>

namespace test {
template <typename Container, typename Allocator>
  struct replace_allocator_impl;

  template <typename Container, typename Allocator>
  using replace_allocator = 
    typename replace_allocator_impl<Container, Allocator>::type;

  template <
    typename K, typename H, typename P, typename A,
    template <typename, typename, typename, typename> class Set,
    typename Allocator
  >
  struct replace_allocator_impl<Set<K, H, P, A>, Allocator>
  {
    using type = Set<
      K, H, P, boost::allocator_rebind_t<Allocator, K> >;
  };

  template <
    typename K, typename H, typename T, typename P, typename A,
    template <typename, typename, typename, typename, typename> class Map,
    typename Allocator
  >
  struct replace_allocator_impl<Map<K, T, H, P, A>, Allocator>
  {
    using type = Map<
      K, T, H, P,
      boost::allocator_rebind_t<Allocator, std::pair<K const, T> > >;
  };
} // namespace test

#endif // !defined(BOOST_UNORDERED_TEST_REPLACE_ALLOCATOR)
