// Copyright 2024 Joaquin M Lopez Munoz.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_UNORDERED_ENABLE_STATS

#ifdef BOOST_UNORDERED_CFOA_TESTS
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_node_map.hpp>
#include <boost/unordered/unordered_node_set.hpp>
#include "../cfoa/helpers.hpp"
#else
#include "../helpers/unordered.hpp"
#endif

#include "../helpers/helpers.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/test.hpp"
#include <boost/assert.hpp>
#include <boost/core/make_span.hpp>
#include <cmath>
#include <cstring>

template <class T> struct unequal_allocator
{
  typedef T value_type;

  unequal_allocator(int n = 0): n_{n} {}
  unequal_allocator(unequal_allocator const&) = default;
  unequal_allocator(unequal_allocator&&) = default;

  template <class U>
  unequal_allocator(unequal_allocator<U> const& x): n_{x.n_} {}

  BOOST_ATTRIBUTE_NODISCARD T* allocate(std::size_t n)
  {
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) noexcept { ::operator delete(p); }

  bool operator==(unequal_allocator const& x) const { return n_ == x.n_; }
  bool operator!=(unequal_allocator const& x) const { return n_ != x.n_; }

  int n_;
};

bool esentially_same(double x, double y)
{
  // Some optimizer-related issues in GCC X86 result in last-bit differences
  // on doubles that should otherwise be identical.

  // https://stackoverflow.com/a/253874/213114

  static constexpr double epsilon = 1.0E-6;
  return fabs(x - y) <= ( (fabs(x) > fabs(y) ? fabs(x) : fabs(y)) * epsilon);
}

bool not_esentially_same(double x, double y)
{
  return !esentially_same(x, y);
}

enum check_stats_contition
{
  stats_empty=0,
  stats_full,
  stats_mostly_full // unsuccesful lookups may result in num_comparisons == 0
};

template <class Stats>
void check_stat(const Stats& s, check_stats_contition cond)
{
  switch (cond) {
  case stats_empty:
    BOOST_TEST(esentially_same(s.average, 0.0));
    BOOST_TEST(esentially_same(s.variance, 0.0));
    BOOST_TEST(esentially_same(s.deviation, 0.0));
    break;
  case stats_full:
    BOOST_TEST_GT(s.average, 0.0);
    if(not_esentially_same(s.variance, 0.0)) {
      BOOST_TEST_GT(s.variance, 0.0);
      BOOST_TEST_GT(s.deviation, 0.0);
    }
    break;
  case stats_mostly_full:
    if(not_esentially_same(s.variance, 0.0)) {
      BOOST_TEST_GT(s.average, 0.0);
      BOOST_TEST_GT(s.variance, 0.0);
      BOOST_TEST_GT(s.deviation, 0.0);
    }
    break;
  default:
    break;
  }
}

template <class Stats> void check_stat(const Stats& s1, const Stats& s2)
{
  BOOST_TEST(esentially_same(s1.average, s2.average));
  BOOST_TEST(esentially_same(s1.variance, s2.variance));
  BOOST_TEST(esentially_same(s1.deviation, s2.deviation));
}

template <class Stats>
void check_insertion_stats(const Stats& s, check_stats_contition cond)
{
  switch (cond) {
  case stats_empty:
    BOOST_TEST_EQ(s.count, 0);
    check_stat(s.probe_length, stats_empty);
    break;
  case stats_full:
    BOOST_TEST_NE(s.count, 0);
    check_stat(s.probe_length, stats_full);
    break;
  default:
    BOOST_ASSERT(false); // insertion can't be mostly full
  }
}

template <class Stats>
void check_insertion_stats(const Stats& s1, const Stats& s2)
{
  BOOST_TEST_EQ(s1.count, s2.count);
  check_stat(s1.probe_length, s2.probe_length);
}

template <class Stats>
void check_lookup_stats(const Stats& s, check_stats_contition cond)
{
  check_stat(s.probe_length, cond == stats_empty? stats_empty : stats_full);
  check_stat(s.num_comparisons, cond);
}

template <class Stats>
void check_lookup_stats(const Stats& s1, const Stats& s2)
{
  BOOST_TEST_EQ(s1.count, s2.count);
  check_stat(s1.probe_length, s2.probe_length);
  check_stat(s1.num_comparisons, s2.num_comparisons);
}

template <class Stats>
void check_container_stats(const Stats& s, check_stats_contition cond)
{
  if(cond == stats_mostly_full) {
    BOOST_ASSERT(false); // mostly full only applies to unsuccessful lookup
  }
  check_insertion_stats(s.insertion, cond);
  check_lookup_stats(s.successful_lookup, cond);
  check_lookup_stats(
    s.unsuccessful_lookup,
    cond == stats_empty? stats_empty : stats_mostly_full);
}

template <class Stats>
void check_container_stats(const Stats& s1, const Stats& s2)
{
  check_insertion_stats(s1.insertion, s2.insertion);
  check_lookup_stats(s1.successful_lookup, s2.successful_lookup);
  check_lookup_stats(s1.unsuccessful_lookup, s2.unsuccessful_lookup);
}

template <class Container> void insert_n(Container& c, std::size_t n)
{
#if defined(BOOST_UNORDERED_CFOA_TESTS)
  using value_type = typename Container::value_type;

  test::reset_sequence();
  test::random_values<Container> l(n, test::sequential);
  std::vector<value_type> v(l.begin(), l.end());
  thread_runner(v, [&c](boost::span<value_type> sp) {
    for (auto const& x : sp) {
      c.insert(x);
    }
  });
#else
  test::reset_sequence();
  test::random_values<Container> l(n, test::sequential);
  c.insert(l.begin(), l.end());
#endif
}

template <class Container> void test_stats()
{
  using allocator_type = typename Container::allocator_type;
  using stats = typename Container::stats;

  Container        c;
  const Container& cc = c;

  // Stats initially empty
  stats s = cc.get_stats(); // using cc -> get_stats() is const
  check_container_stats(s, stats_empty);

  // Stats after insertion
  insert_n(c, 10000);
  s = cc.get_stats();
  check_insertion_stats(s.insertion, stats_full); // insertions happened
  check_lookup_stats(s.successful_lookup, stats_empty); // no duplicate values
  check_lookup_stats(
    s.unsuccessful_lookup, stats_mostly_full); // from insertion

#if !defined(BOOST_UNORDERED_CFOA_TESTS)
  // Inequality due to rehashing
  // May not hold in concurrent containers because of insertion retries
  BOOST_TEST_GT(
    s.insertion.count, s.unsuccessful_lookup.count); 
#endif

  // resets_stats() actually clears stats
  c.reset_stats();
  check_container_stats(cc.get_stats(), stats_empty);

  // Stats after lookup

  test::reset_sequence();

#if defined(BOOST_UNORDERED_CFOA_TESTS)

  using key_type = typename Container::key_type;
  using value_type = typename Container::value_type;

  test::random_values<Container> l2(15000, test::sequential);
  std::vector<value_type>        v2(l2.begin(), l2.end());
  std::atomic<std::size_t>       found{0}, not_found{0};
  thread_runner(v2, [&cc, &found, &not_found](boost::span<value_type> sp) {
    // Half the span looked up elementwise
    auto sp1 = boost::make_span(sp.begin(), sp.size()/2);
    for (auto const& x : sp1) {
      if(cc.contains(test::get_key<Container>(x))) ++found;
      else                                         ++not_found;
    }

    // Second half looked up in bulk
    std::vector<key_type> ksp2;
    for (auto const& x : boost::make_span(
          sp1.end(), static_cast<std::size_t>(sp.end() - sp1.end()))) {
      ksp2.push_back(test::get_key<Container>(x));
    }
    auto visited = cc.visit(
      ksp2.begin(), ksp2.end(), [](const value_type&) {});
    found += visited;
    not_found += ksp2.size() - visited;
  });

#else

  test::random_values<Container> v2(15000, test::sequential);
  std::size_t                    found = 0, not_found = 0;
  for (const auto& x: v2) {
    if (cc.contains(test::get_key<Container>(x))) ++found;
    else                                          ++not_found;
  }

#endif

  // As many [un]successful lookups as recorded externally
  s=cc.get_stats();
  check_lookup_stats(s.successful_lookup, stats_full);
  check_lookup_stats(s.unsuccessful_lookup, stats_mostly_full);
  BOOST_TEST_EQ(s.successful_lookup.count, found);
  BOOST_TEST_EQ(s.unsuccessful_lookup.count, not_found);

  c.reset_stats();
  s = cc.get_stats();
  check_container_stats(s, stats_empty);

  // Move constructor tests

  c.clear();
  insert_n(c, 1000);
  insert_n(c, 1000); // produces successful lookups

  // Move contructor
  // Stats transferred to target and reset in source
  s = cc.get_stats();
  Container c2 = std::move(c);
  check_container_stats(c.get_stats(), stats_empty);
  check_container_stats(c2.get_stats(), s);

  // Move constructor with equal allocator
  // Stats transferred to target and reset in source
  Container c3(std::move(c2), allocator_type());
  check_container_stats(c2.get_stats(), stats_empty);
  check_container_stats(c3.get_stats(), s);

  // Move constructor with unequal allocator
  // Target only has insertions, stats reset in source
  Container c4(std::move(c3), allocator_type(1));
  check_container_stats(c3.get_stats(), stats_empty);
  check_insertion_stats(c4.get_stats().insertion, stats_full);
  check_lookup_stats(c4.get_stats().successful_lookup, stats_empty);
  check_lookup_stats(c4.get_stats().unsuccessful_lookup, stats_empty);

  // Move assignment tests

  // Move assignment with equal allocator
  // Stats transferred to target and reset in source
  Container c5, c6;
  insert_n(c5,1000);
  insert_n(c5,1000); // produces successful lookups
  insert_n(c6,500);
  insert_n(c6,500); // produces successful lookups
  s = c5.get_stats();
  check_container_stats(s, stats_full);
  check_container_stats(c6.get_stats(), stats_full);
  c6 = std::move(c5);
  check_container_stats(c5.get_stats(), stats_empty);
  check_container_stats(c6.get_stats(), s);

  // Move assignment with unequal allocator
  // Target only has insertions (if reset previously), stats reset in source
  Container c7(allocator_type(1));
  insert_n(c7,250);
  insert_n(c7,250); // produces successful lookups
  check_container_stats(c7.get_stats(), stats_full);
  c7.reset_stats();
  c7 = std::move(c6);
  check_container_stats(c6.get_stats(), stats_empty);
  check_insertion_stats(c7.get_stats().insertion, stats_full);
  check_lookup_stats(c7.get_stats().successful_lookup, stats_empty);
  check_lookup_stats(c7.get_stats().unsuccessful_lookup, stats_empty);
}

#if defined(BOOST_UNORDERED_CFOA_TESTS)
template <class Container, class ConcurrentContainer>
void test_stats_concurrent_unordered_interop()
{
  ConcurrentContainer cc1;
  insert_n(cc1,5000);
  insert_n(cc1,5000); // produces successful lookups
  auto s=cc1.get_stats();
  Container c1(std::move(cc1));
  check_container_stats(cc1.get_stats(),stats_empty);
  check_container_stats(c1.get_stats(),s);

  ConcurrentContainer cc2(std::move(c1));
  check_container_stats(c1.get_stats(),stats_empty);
  check_container_stats(cc2.get_stats(),s);
}
#endif

UNORDERED_AUTO_TEST (stats_) {
#if defined(BOOST_UNORDERED_CFOA_TESTS)
  test_stats<
    boost::concurrent_flat_map<
      int, int, boost::hash<int>, std::equal_to<int>,
      unequal_allocator< std::pair< const int, int> >>>();
  test_stats<
    boost::concurrent_node_map<
      int, int, boost::hash<int>, std::equal_to<int>,
      unequal_allocator< std::pair< const int, int> >>>();
  test_stats<
    boost::concurrent_flat_set<
      int, boost::hash<int>, std::equal_to<int>, unequal_allocator<int>>>();
  test_stats<
    boost::concurrent_node_set<
      int, boost::hash<int>, std::equal_to<int>, unequal_allocator<int>>>();
  test_stats_concurrent_unordered_interop<
    boost::unordered_flat_map<int, int>,
    boost::concurrent_flat_map<int, int>>();
  test_stats_concurrent_unordered_interop<
    boost::unordered_node_map<int, int>,
    boost::concurrent_node_map<int, int>>();
  test_stats_concurrent_unordered_interop<
    boost::unordered_flat_set<int>,
    boost::concurrent_flat_set<int>>();
  test_stats_concurrent_unordered_interop<
    boost::unordered_node_set<int>,
    boost::concurrent_node_set<int>>();
#elif defined(BOOST_UNORDERED_FOA_TESTS)
  test_stats<
    boost::unordered_flat_map<
      int, int, boost::hash<int>, std::equal_to<int>,
      unequal_allocator< std::pair< const int, int> >>>();
  test_stats<
    boost::unordered_flat_set<
      int, boost::hash<int>, std::equal_to<int>, unequal_allocator<int>>>();
  test_stats<
    boost::unordered_node_map<
      int, int, boost::hash<int>, std::equal_to<int>,
      unequal_allocator< std::pair< const int, int> >>>();
  test_stats<
    boost::unordered_node_set<
      int, boost::hash<int>, std::equal_to<int>, unequal_allocator<int>>>();
#else
  // Closed-addressing containers do not provide stats
#endif
}

RUN_TESTS()
