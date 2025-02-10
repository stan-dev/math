// Copyright 2021-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"

#include <boost/config.hpp>
#include <boost/container_hash/hash.hpp>

#include <cmath>
#include <functional>

std::size_t total_allocation = 0;
std::size_t num_allocations = 0;

template <typename T> struct A
{
  typedef T value_type;

  static int count;
  int i;

  A() : i(++count) {}

  template <class U> A(const A<U>& a) noexcept : i(a.i) {}

  T* allocate(std::size_t n)
  {
    total_allocation += n * sizeof(T);
    ++num_allocations;
    return (T*)(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t n) noexcept
  {
    total_allocation -= n * sizeof(T);
    ::operator delete(p);
  }

  bool operator==(A const& a) const { return i == a.i; }
  bool operator!=(A const& a) const { return i != a.i; }
};

template <class T> int A<T>::count = 0;

template <class UnorderedContainer>
void bucket_count_constructor(UnorderedContainer*)
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

template <class UnorderedContainer>
void range_bucket_constructor(UnorderedContainer*)
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

template <class UnorderedContainer> void reserve_tests(UnorderedContainer*)
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

    std::size_t expected_bucket_count = static_cast<std::size_t>(
      std::ceil(static_cast<float>(count) / s.max_load_factor()));

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

template <class UnorderedContainer> void rehash_tests(UnorderedContainer*)
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
    std::size_t prev_total_allocation = total_allocation;

    s.rehash(count);
    BOOST_TEST_EQ(num_allocations, prev_allocations);
    BOOST_TEST_EQ(total_allocation, prev_total_allocation);

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
#ifdef BOOST_UNORDERED_FOA_TESTS
    BOOST_TEST_LT(s.bucket_count(), prev_count + count);
    BOOST_TEST_LE(total_allocation,
      (prev_count + count) * sizeof(typename UnorderedContainer::value_type) +
        ((prev_count + count) / 15 + 1) * 16);
#else
    std::size_t const estimated_bucket_group_size =
      3 * sizeof(void*) + sizeof(std::size_t);
    std::size_t const estimated_bucket_groups =
      s.bucket_count() / (sizeof(std::size_t) * 8);

    BOOST_TEST_LT(s.bucket_count(), prev_count + count);
    BOOST_TEST_LE(total_allocation,
      (prev_count + count) * sizeof(void*) +
        estimated_bucket_group_size * estimated_bucket_groups);
#endif
  }

  BOOST_TEST_GT(num_allocations, 0u);
  BOOST_TEST_EQ(total_allocation, 0u);
  num_allocations = 0;
}

UNORDERED_AUTO_TEST (allocator_check) {
  // prove Allocator invariants
  // from cppref:
  // Given:
  // * A, an Allocator type for type T
  // * B, the corresponding Allocator type for some cv-unqualified object type
  //   U (as obtained by rebinding A)
  //
  // Expression:
  // A a(b)
  //
  // Return Type:
  // Constructs `a` such that `B(a)==b` and `A(b)==a`.
  // (Note: This implies that all allocators related by rebind maintain each
  // other's resources, such as memory pools.)
  //
  //
  typedef boost::allocator_rebind<A<int>, float>::type alloc_rebound;
  alloc_rebound b;
  A<int> a(b);
  BOOST_TEST(alloc_rebound(a) == b);
  BOOST_TEST(A<int>(b) == a);
}

#ifdef BOOST_UNORDERED_FOA_TESTS
static boost::unordered_flat_set<int*, boost::hash<int*>, std::equal_to<int*>,
  A<int*> >* test_set;
static boost::unordered_flat_map<int*, int*, boost::hash<int*>,
  std::equal_to<int*>, A<std::pair<int const*, int*> > >* test_map;

static boost::unordered_node_set<int*, boost::hash<int*>, std::equal_to<int*>,
  A<int*> >* test_node_set;

static boost::unordered_node_map<int*, int*, boost::hash<int*>,
  std::equal_to<int*>, A<std::pair<int const*, int*> > >* test_node_map;

// clang-format off
UNORDERED_TEST(bucket_count_constructor,
  ((test_set)(test_map)(test_node_set)(test_node_map)))

UNORDERED_TEST(range_bucket_constructor,
  ((test_set)(test_map)(test_node_set)(test_node_map)))

UNORDERED_TEST(reserve_tests,
  ((test_set)(test_map)(test_node_set)(test_node_map)))

UNORDERED_TEST(rehash_tests,
  ((test_set)(test_map)(test_node_set)(test_node_map)))
// clang-format on
#else
static boost::unordered_set<int, boost::hash<int>, std::equal_to<int>, A<int> >*
  test_set;

static boost::unordered_multiset<int, boost::hash<int>, std::equal_to<int>,
  A<int> >* test_multiset;

static boost::unordered_map<int, int, boost::hash<int>, std::equal_to<int>,
  A<std::pair<int const, int> > >* test_map;

static boost::unordered_multimap<int, int, boost::hash<int>, std::equal_to<int>,
  A<std::pair<int const, int> > >* test_multimap;

// clang-format off
UNORDERED_TEST(bucket_count_constructor,
  ((test_set)(test_map)(test_multiset)(test_multimap)))

UNORDERED_TEST(range_bucket_constructor,
  ((test_set)(test_map)(test_multiset)(test_multimap)))

UNORDERED_TEST(reserve_tests,
  ((test_set)(test_map)(test_multiset)(test_multimap)))

UNORDERED_TEST(rehash_tests,
  ((test_set)(test_map)(test_multiset)(test_multimap)))
// clang-format on
#endif

RUN_TESTS()
