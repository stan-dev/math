// Copyright 2021 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// clang-format off
#include "../helpers/prefix.hpp"
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "../helpers/postfix.hpp"
// clang-format on

#include "../helpers/test.hpp"

#include <boost/config.hpp>
#include <boost/container_hash/hash.hpp>

#include <cmath>
#include <functional>

std::size_t total_allocation = 0;
std::size_t num_allocations = 0;

struct B
{
  B() : i(++count) {}
  static int count;
  int i;
  bool operator==(B const& b) const { return i == b.i; };
  bool operator!=(B const& b) const { return i != b.i; };
};

int B::count = 0;

template <typename T> struct A : B
{
  typedef T value_type;

  A() {}

  template <class U> A(const A<U>&) BOOST_NOEXCEPT {}

  T* allocate(std::size_t n)
  {
    total_allocation += n * sizeof(T);
    ++num_allocations;
    return (T*)std::calloc(n, sizeof(T));
  }

  void deallocate(T* p, std::size_t n) BOOST_NOEXCEPT
  {
    total_allocation -= n * sizeof(T);
    std::free(p);
  }
};

template <class UnorderedContainer> void bucket_count_constructor()
{
  BOOST_TEST_EQ(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);

  {
    std::size_t count = 50000;

    UnorderedContainer s(count);

    BOOST_TEST_GE(total_allocation, count * sizeof(void*));
    BOOST_TEST_GE(s.bucket_count(), count);
  }

  BOOST_TEST_GT(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);
  num_allocations = 0;
}

template <class UnorderedContainer> void range_bucket_constructor()
{
  BOOST_TEST_EQ(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);

  {
    UnorderedContainer s1;

    std::size_t count = 50000;

    UnorderedContainer s2(s1.begin(), s1.end(), count);

    BOOST_TEST_GE(total_allocation, count * sizeof(void*));
    BOOST_TEST_GE(s2.bucket_count(), count);
  }

  BOOST_TEST_GT(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);
  num_allocations = 0;
}

template <class UnorderedContainer> void reserve_tests()
{
  BOOST_TEST_EQ(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);

  {
    UnorderedContainer s;

    // simple math for the test:
    // max_load_factor = max_size / bucket_count, before a rehashing occurs
    //
    // reserve() respects max load factor and its argument implies the max size
    //
    // reserve(count) => bucket_count = ceil(count / max_load_factor)
    // internal policies reshape bucket_count accordingly but guarantee count as
    // a minimum
    //

    std::size_t count = 50000;

    s.max_load_factor(0.37f);
    s.reserve(count);

    std::size_t expected_bucket_count =
      static_cast<std::size_t>(std::ceil(static_cast<float>(count) / 0.37f));

    BOOST_TEST_GE(total_allocation, expected_bucket_count * sizeof(void*));
    BOOST_TEST_GE(s.bucket_count(), expected_bucket_count);

    std::size_t prev_allocations = num_allocations;
    s.reserve(count);
    BOOST_TEST_EQ(num_allocations, prev_allocations);
  }

  BOOST_TEST_GT(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);
  num_allocations = 0;
}

template <class UnorderedContainer> void rehash_tests()
{
  BOOST_TEST_EQ(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);

  {
    UnorderedContainer s;

    std::size_t count = 1000;
    s.rehash(count);

    // test that an initial allocation occurs
    //
    BOOST_TEST_GE(total_allocation, count * sizeof(void*));
    BOOST_TEST_GE(s.bucket_count(), count);

    // prove idempotence, that rehashing with the exact same bucket count causes
    // no reallocations
    //
    std::size_t prev_allocations = num_allocations;

    s.rehash(count);
    BOOST_TEST_EQ(num_allocations, prev_allocations);

    // prove that when we rehash, exceeding the current bucket count, that we
    // properly deallocate the current bucket array and then reallocate the
    // larger one
    //
    std::size_t prev_count = s.bucket_count();

    count = s.bucket_count() + 2;
    s.rehash(count);

    BOOST_TEST_GT(num_allocations, prev_allocations);
    BOOST_TEST_GE(total_allocation, count * sizeof(void*));
    BOOST_TEST_GE(s.bucket_count(), count);

    // concurrent memory usage here should be less than the sum of the memory
    // required for the previous bucket array and our current one
    // note, the test is vulnerable to cases where the next calculated bucket
    // count can exceed `prev_count + count`
    //
    BOOST_TEST_LT(s.bucket_count(), prev_count + count);
    BOOST_TEST_LT(total_allocation, (prev_count + count) * sizeof(void*));
  }

  BOOST_TEST_GT(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);
  num_allocations = 0;
}

UNORDERED_AUTO_TEST (unordered_set_reserve) {
  typedef boost::unordered_set<int, boost::hash<int>, std::equal_to<int>,
    A<int> >
    unordered_set;

  typedef boost::unordered_multiset<int, boost::hash<int>, std::equal_to<int>,
    A<int> >
    unordered_multiset;

  typedef boost::unordered_map<int, int, boost::hash<int>, std::equal_to<int>,
    A<std::pair<int const, int> > >
    unordered_map;

  typedef boost::unordered_multimap<int, int, boost::hash<int>,
    std::equal_to<int>, A<std::pair<int const, int> > >
    unordered_multimap;

  bucket_count_constructor<unordered_set>();
  bucket_count_constructor<unordered_map>();
  bucket_count_constructor<unordered_multiset>();
  bucket_count_constructor<unordered_multimap>();

  range_bucket_constructor<unordered_set>();
  range_bucket_constructor<unordered_map>();
  range_bucket_constructor<unordered_multiset>();
  range_bucket_constructor<unordered_multimap>();

  reserve_tests<unordered_set>();
  reserve_tests<unordered_map>();
  reserve_tests<unordered_multiset>();
  reserve_tests<unordered_multimap>();

  rehash_tests<unordered_set>();
  rehash_tests<unordered_map>();
  rehash_tests<unordered_multiset>();
  rehash_tests<unordered_multimap>();
}

RUN_TESTS()
