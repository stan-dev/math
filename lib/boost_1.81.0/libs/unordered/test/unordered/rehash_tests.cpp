
// Copyright 2006-2009 Daniel James.
// Copyright 2022 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/metafunctions.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/test.hpp"
#include "../helpers/tracker.hpp"
#include "../objects/test.hpp"

namespace rehash_tests {

  test::seed_t initialize_seed(2974);

  static int count_allocations;
  template <class T> struct monotonic_allocator
  {
    typedef T value_type;
    monotonic_allocator() {}
    monotonic_allocator(monotonic_allocator const&) {}

    template <class U> monotonic_allocator(monotonic_allocator<U> const&) {}

    friend bool operator==(
      monotonic_allocator const&, monotonic_allocator const&)
    {
      return true;
    }

    friend bool operator!=(
      monotonic_allocator const&, monotonic_allocator const&)
    {
      return false;
    }

    T* allocate(std::size_t n)
    {
      ++count_allocations;
      return static_cast<T*>(::operator new(sizeof(T) * n));
    }

    void deallocate(T* p, std::size_t) { ::operator delete(p); }
  };

  void reset_counts() { count_allocations = 0; }

  template <class X> bool postcondition(X const& x, typename X::size_type n)
  {
    return static_cast<double>(x.bucket_count()) >=
             static_cast<double>(x.size()) / x.max_load_factor() &&
           x.bucket_count() >= n;
  }

  template <class X> void rehash_empty_test1(X*)
  {
    X x;

    x.rehash(10000);
    BOOST_TEST(postcondition(x, 10000));

    x.rehash(0);
    BOOST_TEST(postcondition(x, 0));

    x.rehash(10000000);
    BOOST_TEST(postcondition(x, 10000000));
  }

  template <class X>
  void rehash_empty_test2(X*, test::random_generator generator)
  {
    test::random_values<X> v(1000, generator);
    test::ordered<X> tracker;

    X x;

    x.rehash(10000);
    BOOST_TEST(postcondition(x, 10000));

    tracker.insert_range(v.begin(), v.end());
    x.insert(v.begin(), v.end());
    tracker.compare(x);

    BOOST_TEST(postcondition(x, 10000));

    x.rehash(10000000);
    tracker.compare(x);
    BOOST_TEST(postcondition(x, 10000000));
  }

  template <class X>
  void rehash_empty_test3(X*, test::random_generator generator)
  {
    test::random_values<X> v(1000, generator);
    test::ordered<X> tracker;

    X x;

    x.rehash(0);
    BOOST_TEST(postcondition(x, 0));

    tracker.insert_range(v.begin(), v.end());
    x.insert(v.begin(), v.end());
    tracker.compare(x);

    BOOST_TEST(postcondition(x, 0));
  }

  template <class X> void rehash_empty_tracking(X*, test::random_generator)
  {
    // valid for all load factors
    float const max_load_factors[] = {
      0.5f, 1.0f, 1e6f, std::numeric_limits<float>::infinity()};

    std::size_t const max_load_factors_len =
      sizeof(max_load_factors) / sizeof(*max_load_factors);

    for (std::size_t i = 0; i < max_load_factors_len; ++i) {
      X x;
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);

      x.max_load_factor(max_load_factors[i]);

      {
        BOOST_TEST_EQ(x.bucket_count(), 0u);

        x.rehash(0);
        BOOST_TEST_EQ(x.bucket_count(), 0u);
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }

      {
        BOOST_TEST_EQ(x.bucket_count(), 0u);

        x.rehash(1000);
        BOOST_TEST_GE(x.bucket_count(), 1000u);
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        x.rehash(0);
        BOOST_TEST_EQ(x.bucket_count(), 0u);
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }

      {
        BOOST_TEST_EQ(x.bucket_count(), 0u);

        x.rehash(1000);
        BOOST_TEST_GE(x.bucket_count(), 1000u);
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        x.rehash(10);
        BOOST_TEST_GE(x.bucket_count(), 10u);
        BOOST_TEST_LT(x.bucket_count(), 1000u);
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);
      }

      {
        BOOST_TEST_GT(x.bucket_count(), 0u);
        BOOST_TEST_LT(x.bucket_count(), 1000u);

        x.rehash(1000);
        BOOST_TEST_GE(x.bucket_count(), 1000u);
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        x.rehash(0);
        BOOST_TEST_EQ(x.bucket_count(), 0u);
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }
    }

    for (std::size_t i = 0; i < max_load_factors_len; ++i) {
      typedef typename X::size_type size_type;

      X x;
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);

#ifdef BOOST_UNORDERED_FOA_TESTS
      float const mlf = boost::unordered::detail::foa::mlf;
      x.max_load_factor(max_load_factors[i]);
#else
      float const mlf = max_load_factors[i];
      x.max_load_factor(mlf);
#endif

      {
        BOOST_TEST_EQ(x.bucket_count(), 0u);

        x.reserve(0);
        BOOST_TEST_EQ(x.bucket_count(), 0u);
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }

      {
        BOOST_TEST_EQ(x.bucket_count(), 0u);

        x.reserve(1000);
        BOOST_TEST_GE(
          x.bucket_count(), static_cast<size_type>(std::ceil(1000 / mlf)));
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        x.reserve(0);
        BOOST_TEST_EQ(x.bucket_count(), 0u);
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }

      {
        BOOST_TEST_EQ(x.bucket_count(), 0u);

        x.reserve(1000);
        BOOST_TEST_GE(
          x.bucket_count(), static_cast<size_type>(std::ceil(1000 / mlf)));
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        x.reserve(10);
        BOOST_TEST_GE(
          x.bucket_count(), static_cast<size_type>(std::ceil(10 / mlf)));
        BOOST_TEST_LT(x.bucket_count(), 1000u);
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);
      }

      {
        BOOST_TEST_GT(x.bucket_count(), 0u);
        BOOST_TEST_LT(x.bucket_count(), 1000u);

        x.reserve(1000);
        BOOST_TEST_GE(
          x.bucket_count(), static_cast<size_type>(std::ceil(1000 / mlf)));
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        x.reserve(0);
        BOOST_TEST_EQ(x.bucket_count(), 0u);
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }
    }
  }

  template <class X>
  void rehash_nonempty_tracking(X*, test::random_generator generator)
  {
    test::random_values<X> const v(1000, generator);

    typedef typename X::size_type size_type;

    float const max_load_factors[] = {0.5f, 1.0f, 1e2f};

    size_type const max_load_factors_len =
      sizeof(max_load_factors) / sizeof(*max_load_factors);

    for (size_type i = 0; i < max_load_factors_len; ++i) {
      float const mlf = max_load_factors[i];

      X x(v.begin(), v.end());
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

      x.max_load_factor(mlf);

      size_type bucket_count = x.bucket_count();

      {
        BOOST_TEST_GT(x.bucket_count(), 0u);

        x.rehash(0);
        BOOST_TEST_GE(x.bucket_count(),
          static_cast<size_type>(
            std::floor(static_cast<float>(x.size()) / x.max_load_factor())));
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        bucket_count = x.bucket_count();
      }

      {
        BOOST_TEST_GT(bucket_count, 0u);

        x.rehash(2 * x.bucket_count());
        BOOST_TEST_GT(x.bucket_count(), bucket_count);

        bucket_count = x.bucket_count();
      }

      {
        float const old_mlf = x.max_load_factor();

        BOOST_TEST_GT(bucket_count, 0u);

        x.rehash(bucket_count / 4);
        BOOST_TEST_LT(x.bucket_count(), bucket_count);

        x.max_load_factor(std::numeric_limits<float>::infinity());
        x.rehash(0);
        BOOST_TEST_GT(x.bucket_count(), 0u);

        x.max_load_factor(old_mlf);
      }

      {
        std::size_t const max_load =
          static_cast<std::size_t>(static_cast<double>(x.max_load_factor()) *
                                   static_cast<double>(x.bucket_count()));

        while (x.size() < max_load) {
          test::random_values<X> const t(max_load, generator);
          typename test::random_values<X>::const_iterator pos = t.begin();
          typename test::random_values<X>::const_iterator end = t.end();
          for (; pos != end; ++pos) {
            x.insert(*pos);
            if (x.size() == max_load) {
              break;
            }
          }
        }

        while (x.size() > max_load) {
          x.erase(x.begin());
        }

        BOOST_TEST_EQ(x.size(), max_load);

        bucket_count = x.bucket_count();
        x.rehash(x.bucket_count());
        BOOST_TEST_EQ(x.bucket_count(), bucket_count);
      }
    }

    for (size_type i = 0; i < max_load_factors_len; ++i) {
      X x(v.begin(), v.end());
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

      float const mlf = max_load_factors[i];
      x.max_load_factor(mlf);

      size_type bucket_count = x.bucket_count();

      {
        BOOST_TEST_GT(x.bucket_count(), 0u);

        x.reserve(0);
        BOOST_TEST_GE(x.bucket_count(),
          static_cast<size_type>(
            std::floor(static_cast<float>(x.size()) / x.max_load_factor())));
        BOOST_TEST_GT(test::detail::tracker.count_allocations, 0u);

        bucket_count = x.bucket_count();
      }

      {
        BOOST_TEST_GT(x.bucket_count(), 0u);

        x.reserve(
          2 *
          (static_cast<size_type>(
            std::floor(static_cast<float>(x.size()) / x.max_load_factor()) +
            std::floor(static_cast<float>(x.size()) * x.max_load_factor()))));

        BOOST_TEST_GT(x.bucket_count(), bucket_count);

        bucket_count = x.bucket_count();
        BOOST_TEST_GT(bucket_count, 1u);
      }

      {
        float const old_mlf = x.max_load_factor();

        BOOST_TEST_GT(bucket_count, 4u);

        x.reserve(bucket_count / 4);
        BOOST_TEST_LT(x.bucket_count(), bucket_count);

        x.max_load_factor(std::numeric_limits<float>::infinity());
        x.reserve(0);
        BOOST_TEST_GT(x.bucket_count(), 0u);

        x.max_load_factor(old_mlf);
      }

      {
        std::size_t const max_load =
          static_cast<std::size_t>(static_cast<double>(x.max_load_factor()) *
                                   static_cast<double>(x.bucket_count()));

        while (x.size() < max_load) {
          test::random_values<X> const t(max_load, generator);
          typename test::random_values<X>::const_iterator pos = t.begin();
          typename test::random_values<X>::const_iterator end = t.end();
          for (; pos != end; ++pos) {
            x.insert(*pos);
            if (x.size() == max_load) {
              break;
            }
          }
        }

        while (x.size() > max_load) {
          x.erase(x.begin());
        }

        BOOST_TEST_EQ(x.size(), max_load);

        bucket_count = x.bucket_count();
        x.reserve(x.size());
        BOOST_TEST_EQ(x.bucket_count(), bucket_count);
      }
    }
  }

  template <class X> void rehash_stability(X*, test::random_generator generator)
  {
    reset_counts();

    typedef typename X::size_type size_type;

    size_type bucket_count = 100;
    X x(bucket_count);

    size_type num_elems = x.bucket_count() - 1;

    test::random_values<X> v(num_elems, generator);
    test::ordered<X> tracker;
    tracker.insert_range(v.begin(), v.end());

    typename test::random_values<X>::iterator pos = v.begin();
    for (size_type i = 0; i < num_elems; ++i) {
      x.insert(*pos);
      ++pos;
    }

    int const old_count = count_allocations;
    x.rehash(0);

    BOOST_TEST_EQ(count_allocations, old_count);
    tracker.compare(x);
  }

  template <class X> void rehash_test1(X*, test::random_generator generator)
  {
    test::random_values<X> v(1000, generator);
    test::ordered<X> tracker;
    tracker.insert_range(v.begin(), v.end());
    X x(v.begin(), v.end());

    x.rehash(0);
    BOOST_TEST(postcondition(x, 0));
    tracker.compare(x);

    x.max_load_factor(0.25);
    x.rehash(0);
    BOOST_TEST(postcondition(x, 0));
    tracker.compare(x);

    x.max_load_factor(50.0);
    x.rehash(0);
    BOOST_TEST(postcondition(x, 0));
    tracker.compare(x);

    x.rehash(1000);
    BOOST_TEST(postcondition(x, 1000));
    tracker.compare(x);
  }

  template <class X> void reserve_empty_test1(X*)
  {
    X x;

    x.reserve(10000);
    BOOST_TEST(x.bucket_count() >= 10000);

    x.reserve(0);

    x.reserve(10000000);
    BOOST_TEST(x.bucket_count() >= 10000000);
  }

  template <class X> void reserve_empty_test2(X*)
  {
    X x;
    x.max_load_factor(0.25);

#ifdef BOOST_UNORDERED_FOA_TESTS
    x.reserve(10000);
    BOOST_TEST(x.bucket_count() >= 10000);

    x.reserve(0);

    x.reserve(10000000);
    BOOST_TEST(x.bucket_count() >= 10000000);
#else
    x.reserve(10000);
    BOOST_TEST(x.bucket_count() >= 40000);

    x.reserve(0);

    x.reserve(10000000);
    BOOST_TEST(x.bucket_count() >= 40000000);
#endif
  }

  template <class X> void reserve_test1(X*, test::random_generator generator)
  {
    for (int random_mlf = 0; random_mlf < 2; ++random_mlf) {
      for (std::size_t i = 1; i < 2000; i += i < 50 ? 1 : 13) {
        test::random_values<X> v(i, generator);

        test::ordered<X> tracker;
        tracker.insert_range(v.begin(), v.end());

        X x;
        x.max_load_factor(
          random_mlf ? static_cast<float>(std::rand() % 1000) / 500.0f + 0.5f
                     : 1.0f);
        x.reserve(test::has_unique_keys<X>::value ? i : v.size());

        // Insert an element before the range insert, otherwise there are
        // no iterators to invalidate in the range insert, and it can
        // rehash.
        typename test::random_values<X>::iterator it = v.begin();
        x.insert(*it);
        ++it;

        std::size_t bucket_count = x.bucket_count();
        x.insert(it, v.end());
        BOOST_TEST(bucket_count == x.bucket_count());
        tracker.compare(x);
      }
    }
  }

  template <class X> void reserve_test2(X*, test::random_generator generator)
  {
    for (int random_mlf = 0; random_mlf < 2; ++random_mlf) {
      for (std::size_t i = 0; i < 2000; i += i < 50 ? 1 : 13) {
        test::random_values<X> v(i, generator);

        test::ordered<X> tracker;
        tracker.insert_range(v.begin(), v.end());

        X x;
        x.max_load_factor(
          random_mlf ? static_cast<float>(std::rand() % 1000) / 500.0f + 0.5f
                     : 1.0f);

        x.reserve(test::has_unique_keys<X>::value ? i : v.size());

        std::size_t bucket_count = x.bucket_count();
        for (typename test::random_values<X>::iterator it = v.begin();
             it != v.end(); ++it) {
          x.insert(*it);
        }

        BOOST_TEST(bucket_count == x.bucket_count());
        tracker.compare(x);
      }
    }
  }

  using test::default_generator;
  using test::generate_collisions;
  using test::limited_range;

#ifdef BOOST_UNORDERED_FOA_TESTS
  boost::unordered_flat_set<int>* int_set_ptr;
  boost::unordered_flat_map<test::movable, test::movable, test::hash,
    test::equal_to, test::allocator2<test::movable> >* test_map_ptr;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_set_tracking;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    test::allocator1<std::pair<test::object const, test::object> > >*
    test_map_tracking;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    monotonic_allocator<test::object> >* test_set_monotonic;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    monotonic_allocator<std::pair<test::object const, test::object> > >*
    test_map_monotonic;

  UNORDERED_TEST(rehash_empty_test1, ((int_set_ptr)(test_map_ptr)))
  UNORDERED_TEST(rehash_empty_test2,
    ((int_set_ptr)(test_map_ptr))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(rehash_empty_test3,
    ((int_set_ptr)(test_map_ptr))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(
    rehash_test1, ((int_set_ptr)(test_map_ptr))(
                    (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(reserve_empty_test1, ((int_set_ptr)(test_map_ptr)))
  UNORDERED_TEST(reserve_empty_test2, ((int_set_ptr)(test_map_ptr)))
  UNORDERED_TEST(
    reserve_test1, ((int_set_ptr)(test_map_ptr))(
                     (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(
    reserve_test2, ((int_set_ptr)(test_map_ptr))(
                     (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(rehash_empty_tracking,
    ((test_set_tracking)(test_map_tracking))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(
    rehash_nonempty_tracking, ((test_set_tracking)(test_map_tracking))(
                                (default_generator)(limited_range)))
  UNORDERED_TEST(rehash_stability, ((test_set_monotonic)(test_map_monotonic))(
                                     (default_generator)(limited_range)))
#else
  boost::unordered_set<int>* int_set_ptr;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::allocator2<test::object> >* test_multiset_ptr;
  boost::unordered_map<test::movable, test::movable, test::hash, test::equal_to,
    test::allocator2<test::movable> >* test_map_ptr;
  boost::unordered_multimap<int, int>* int_multimap_ptr;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_set_tracking;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_multiset_tracking;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    test::allocator1<std::pair<test::object const, test::object> > >*
    test_map_tracking;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to,
    test::allocator1<std::pair<test::object const, test::object> > >*
    test_multimap_tracking;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    monotonic_allocator<test::object> >* test_set_monotonic;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    monotonic_allocator<test::object> >* test_multiset_monotonic;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    monotonic_allocator<std::pair<test::object const, test::object> > >*
    test_map_monotonic;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to,
    monotonic_allocator<std::pair<test::object const, test::object> > >*
    test_multimap_monotonic;

  UNORDERED_TEST(rehash_empty_test1,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr)))
  UNORDERED_TEST(rehash_empty_test2,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(rehash_empty_test3,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(rehash_test1,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(reserve_empty_test1,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr)))
  UNORDERED_TEST(reserve_empty_test2,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr)))
  UNORDERED_TEST(reserve_test1,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(reserve_test2,
    ((int_set_ptr)(test_multiset_ptr)(test_map_ptr)(int_multimap_ptr))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(rehash_empty_tracking,
    ((test_set_tracking)(test_multiset_tracking)(test_map_tracking)(test_multimap_tracking))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(rehash_nonempty_tracking,
    ((test_set_tracking)(test_multiset_tracking)(test_map_tracking)(test_multimap_tracking))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(rehash_stability,
    ((test_set_monotonic)(test_multiset_monotonic)(test_map_monotonic)(test_multimap_monotonic))(
      (default_generator)(limited_range)))
#endif
} // namespace rehash_tests

RUN_TESTS()
