// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

test::seed_t initialize_seed(4122023);

using test::default_generator;
using test::limited_range;
using test::sequential;

template <class T> struct soccc_allocator
{
  int x_ = -1;

  using value_type = T;

  soccc_allocator() = default;
  soccc_allocator(soccc_allocator const&) = default;
  soccc_allocator(soccc_allocator&&) = default;

  soccc_allocator(int const x) : x_{x} {}

  template <class U> soccc_allocator(soccc_allocator<U> const& rhs) : x_{rhs.x_}
  {
  }

  T* allocate(std::size_t n)
  {
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) { ::operator delete(p); }

  soccc_allocator select_on_container_copy_construction() const
  {
    return {x_ + 1};
  }

  bool operator==(soccc_allocator const& rhs) const { return x_ == rhs.x_; }
  bool operator!=(soccc_allocator const& rhs) const { return x_ != rhs.x_; }
};

using hasher = stateful_hash;
using key_equal = stateful_key_equal;
using allocator_type = stateful_allocator<std::pair<raii const, raii> >;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, allocator_type>;

using map_value_type = typename map_type::value_type;

UNORDERED_AUTO_TEST (default_constructor) {
  boost::unordered::concurrent_flat_map<raii, raii> x;
  BOOST_TEST(x.empty());
  BOOST_TEST_EQ(x.size(), 0u);
}

UNORDERED_AUTO_TEST (bucket_count_with_hasher_key_equal_and_allocator) {
  raii::reset_counts();
  {
    map_type x(0);

    BOOST_TEST(x.empty());
    BOOST_TEST_EQ(x.size(), 0u);
    BOOST_TEST_EQ(x.hash_function(), hasher());
    BOOST_TEST_EQ(x.key_eq(), key_equal());
  }

  {
    map_type x(0, hasher(1));

    BOOST_TEST(x.empty());
    BOOST_TEST_EQ(x.size(), 0u);
    BOOST_TEST_EQ(x.hash_function(), hasher(1));
    BOOST_TEST_EQ(x.key_eq(), key_equal());
  }

  {
    map_type x(0, hasher(1), key_equal(2));

    BOOST_TEST(x.empty());
    BOOST_TEST_EQ(x.size(), 0u);
    BOOST_TEST_EQ(x.hash_function(), hasher(1));
    BOOST_TEST_EQ(x.key_eq(), key_equal(2));
  }

  {
    map_type x(0, hasher(1), key_equal(2), allocator_type{});

    BOOST_TEST(x.empty());
    BOOST_TEST_EQ(x.size(), 0u);
    BOOST_TEST_EQ(x.hash_function(), hasher(1));
    BOOST_TEST_EQ(x.key_eq(), key_equal(2));
    BOOST_TEST(x.get_allocator() == allocator_type{});
  }
}

UNORDERED_AUTO_TEST (soccc) {
  raii::reset_counts();

  boost::unordered::concurrent_flat_map<raii, raii, hasher, key_equal,
    soccc_allocator<std::pair<raii const, raii> > >
    x;

  boost::unordered::concurrent_flat_map<raii, raii, hasher, key_equal,
    soccc_allocator<std::pair<raii const, raii> > >
    y(x);

  BOOST_TEST_EQ(y.hash_function(), x.hash_function());
  BOOST_TEST_EQ(y.key_eq(), x.key_eq());
  BOOST_TEST(y.get_allocator() != x.get_allocator());
}

namespace {
  template <class G> void from_iterator_range(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());
    raii::reset_counts();

    {
      map_type x(values.begin(), values.end());

      test_matches_reference(x, reference_map);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      map_type x(values.begin(), values.end(), 0);

      test_matches_reference(x, reference_map);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      map_type x(values.begin(), values.end(), 0, hasher(1));

      test_matches_reference(x, reference_map);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2));

      test_matches_reference(x, reference_map);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{});

      test_matches_reference(x, reference_map);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    check_raii_counts();
  }

  template <class G> void copy_constructor(G gen, test::random_generator rg)
  {
    {
      map_type x(0, hasher(1), key_equal(2), allocator_type{});
      map_type y(x);

      BOOST_TEST_EQ(y.size(), x.size());
      BOOST_TEST_EQ(y.hash_function(), x.hash_function());
      BOOST_TEST_EQ(y.key_eq(), x.key_eq());
      BOOST_TEST(y.get_allocator() == x.get_allocator());
    }

    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    {
      map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{});

      thread_runner(
        values, [&x, &reference_map](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;
          map_type y(x);

          test_matches_reference(x, reference_map);
          test_matches_reference(y, reference_map);
          BOOST_TEST_EQ(y.size(), x.size());
          BOOST_TEST_EQ(y.hash_function(), x.hash_function());
          BOOST_TEST_EQ(y.key_eq(), x.key_eq());
          BOOST_TEST(y.get_allocator() == x.get_allocator());
        });
    }

    check_raii_counts();

    raii::reset_counts();

    {
      allocator_type a;

      map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2), a);

      thread_runner(
        values, [&x, &reference_map, a](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;
          map_type y(x, a);

          test_matches_reference(x, reference_map);
          test_matches_reference(y, reference_map);
          BOOST_TEST_EQ(y.size(), x.size());
          BOOST_TEST_EQ(y.hash_function(), x.hash_function());
          BOOST_TEST_EQ(y.key_eq(), x.key_eq());
          BOOST_TEST(y.get_allocator() == x.get_allocator());
        });
    }

    check_raii_counts();
  }

  template <class G>
  void copy_constructor_with_insertion(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());
    raii::reset_counts();

    std::mutex m;
    std::condition_variable cv;
    bool ready = false;

    {
      map_type x(0, hasher(1), key_equal(2), allocator_type{});

      auto f = [&x, &values, &m, &cv, &ready] {
        {
          std::lock_guard<std::mutex> guard(m);
          ready = true;
        }
        cv.notify_all();

        for (auto const& val : values) {
          x.insert(val);
        }
      };

      std::thread t1(f);
      std::thread t2(f);

      thread_runner(
        values, [&x, &reference_map, &values, rg, &m, &cv, &ready](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [&] { return ready; });
          }

          map_type y(x);

          BOOST_TEST_LE(y.size(), values.size());
          BOOST_TEST_EQ(y.hash_function(), x.hash_function());
          BOOST_TEST_EQ(y.key_eq(), x.key_eq());
          BOOST_TEST(y.get_allocator() == x.get_allocator());

          x.visit_all([&reference_map, rg](
                        typename map_type::value_type const& val) {
            BOOST_TEST(reference_map.contains(val.first));
            if (rg == sequential) {
              BOOST_TEST_EQ(val.second, reference_map.find(val.first)->second);
            }
          });
        });

      t1.join();
      t2.join();
    }

    check_raii_counts();
  }

  template <class G> void move_constructor(G gen, test::random_generator rg)
  {
    {
      map_type x(0, hasher(1), key_equal(2), allocator_type{});
      auto const old_size = x.size();

      map_type y(std::move(x));

      BOOST_TEST_EQ(y.size(), old_size);
      BOOST_TEST_EQ(y.hash_function(), hasher(1));
      BOOST_TEST_EQ(y.key_eq(), key_equal(2));

      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());

      BOOST_TEST(y.get_allocator() == x.get_allocator());
    }

    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    {
      map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{});

      std::atomic_uint num_transfers{0};

      auto const old_mc = +raii::move_constructor;

      thread_runner(
        values, [&x, &reference_map, &num_transfers](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          auto const old_size = x.size();
          map_type y(std::move(x));

          if (!y.empty()) {
            ++num_transfers;

            test_matches_reference(y, reference_map);
            BOOST_TEST_EQ(y.size(), old_size);
            BOOST_TEST_EQ(y.hash_function(), hasher(1));
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.size(), 0u);
            BOOST_TEST_EQ(y.hash_function(), hasher());
            BOOST_TEST_EQ(y.key_eq(), key_equal());
          }

          BOOST_TEST_EQ(x.size(), 0u);
          BOOST_TEST_EQ(x.hash_function(), hasher());
          BOOST_TEST_EQ(x.key_eq(), key_equal());

          BOOST_TEST(y.get_allocator() == x.get_allocator());
        });

      BOOST_TEST_EQ(num_transfers, 1u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
    }

    check_raii_counts();

    // allocator-aware move constructor, unequal allocators
    raii::reset_counts();

    {
      map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{1});

      std::atomic_uint num_transfers{0};

      auto const old_mc = +raii::move_constructor;
      auto const old_size = x.size();

      thread_runner(
        values, [&x, &reference_map, &num_transfers, old_size](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          auto a = allocator_type{2};
          BOOST_TEST(a != x.get_allocator());

          map_type y(std::move(x), a);

          if (!y.empty()) {
            ++num_transfers;

            test_matches_reference(y, reference_map);
            BOOST_TEST_EQ(y.size(), old_size);
            BOOST_TEST_EQ(y.hash_function(), hasher(1));
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.size(), 0u);
            BOOST_TEST_EQ(y.hash_function(), hasher());
            BOOST_TEST_EQ(y.key_eq(), key_equal());
          }

          BOOST_TEST_EQ(x.size(), 0u);
          BOOST_TEST_EQ(x.hash_function(), hasher());
          BOOST_TEST_EQ(x.key_eq(), key_equal());

          BOOST_TEST(y.get_allocator() != x.get_allocator());
          BOOST_TEST(y.get_allocator() == a);
        });

      BOOST_TEST_EQ(num_transfers, 1u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc + (2 * old_size));
    }

    check_raii_counts();

    // allocator-aware move constructor, equal allocators
    raii::reset_counts();

    {
      map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{1});

      std::atomic_uint num_transfers{0};

      auto const old_mc = +raii::move_constructor;
      auto const old_size = x.size();

      thread_runner(
        values, [&x, &reference_map, &num_transfers, old_size](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          auto a = allocator_type{1};
          BOOST_TEST(a == x.get_allocator());

          map_type y(std::move(x), a);

          if (!y.empty()) {
            ++num_transfers;

            test_matches_reference(y, reference_map);
            BOOST_TEST_EQ(y.size(), old_size);
            BOOST_TEST_EQ(y.hash_function(), hasher(1));
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.size(), 0u);
            BOOST_TEST_EQ(y.hash_function(), hasher());
            BOOST_TEST_EQ(y.key_eq(), key_equal());
          }

          BOOST_TEST_EQ(x.size(), 0u);
          BOOST_TEST_EQ(x.hash_function(), hasher());
          BOOST_TEST_EQ(x.key_eq(), key_equal());

          BOOST_TEST(y.get_allocator() == x.get_allocator());
          BOOST_TEST(y.get_allocator() == a);
        });

      BOOST_TEST_EQ(num_transfers, 1u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
    }

    check_raii_counts();
  }

  template <class G>
  void move_constructor_with_insertion(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    std::mutex m;
    std::condition_variable cv;
    bool ready = false;

    {
      map_type x(0, hasher(1), key_equal(2), allocator_type{});

      std::atomic_uint num_transfers{0};

      std::thread t1([&x, &values] {
        for (auto const& val : values) {
          x.insert(val);
        }
      });

      std::thread t2([&x, &m, &cv, &ready] {
        while (x.empty()) {
          std::this_thread::yield();
        }

        {
          std::lock_guard<std::mutex> guard(m);
          ready = true;
        }
        cv.notify_all();
      });

      thread_runner(
        values, [&x, &reference_map, &num_transfers, rg, &m, &ready, &cv](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [&] { return ready; });
          }

          map_type y(std::move(x));

          if (!y.empty()) {
            ++num_transfers;
            y.cvisit_all([&reference_map, rg](map_value_type const& val) {
              BOOST_TEST(reference_map.contains(val.first));
              if (rg == sequential) {
                BOOST_TEST_EQ(
                  val.second, reference_map.find(val.first)->second);
              }
            });
          }
        });

      t1.join();
      t2.join();

      BOOST_TEST_GE(num_transfers, 1u);
    }

    check_raii_counts();
  }

  template <class G>
  void iterator_range_with_allocator(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    {
      allocator_type a;
      map_type x(values.begin(), values.end(), a);

      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }

      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());

      BOOST_TEST(x.get_allocator() == a);

      test_fuzzy_matches_reference(x, reference_map, rg);
    }

    check_raii_counts();
  }

  UNORDERED_AUTO_TEST (explicit_allocator) {
    raii::reset_counts();

    {
      allocator_type a;
      map_type x(a);

      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());

      BOOST_TEST(x.get_allocator() == a);
    }
  }

  UNORDERED_AUTO_TEST (initializer_list_with_all_params) {
    // hard-code 11 unique values
    std::initializer_list<map_value_type> ilist{
      map_value_type{raii{0}, raii{0}},
      map_value_type{raii{1}, raii{1}},
      map_value_type{raii{2}, raii{2}},
      map_value_type{raii{3}, raii{3}},
      map_value_type{raii{4}, raii{4}},
      map_value_type{raii{5}, raii{5}},
      map_value_type{raii{6}, raii{6}},
      map_value_type{raii{6}, raii{6}},
      map_value_type{raii{7}, raii{7}},
      map_value_type{raii{8}, raii{8}},
      map_value_type{raii{9}, raii{9}},
      map_value_type{raii{10}, raii{10}},
      map_value_type{raii{9}, raii{9}},
      map_value_type{raii{8}, raii{8}},
      map_value_type{raii{7}, raii{7}},
      map_value_type{raii{6}, raii{6}},
      map_value_type{raii{5}, raii{5}},
      map_value_type{raii{4}, raii{4}},
      map_value_type{raii{3}, raii{3}},
      map_value_type{raii{2}, raii{2}},
      map_value_type{raii{1}, raii{1}},
      map_value_type{raii{0}, raii{0}},
    };

    {
      raii::reset_counts();

      map_type x(ilist, 0, hasher(1), key_equal(2), allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 2 * ilist.size());
      BOOST_TEST_EQ(raii::move_constructor, 2 * 11u);
    }
    check_raii_counts();

    {
      raii::reset_counts();

      map_type x(ilist, allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 2 * ilist.size());
      BOOST_TEST_EQ(raii::move_constructor, 2 * 11u);
    }
    check_raii_counts();

    {
      raii::reset_counts();

      map_type x(ilist, 0, allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 2 * ilist.size());
      BOOST_TEST_EQ(raii::move_constructor, 2 * 11u);
    }
    check_raii_counts();

    {
      raii::reset_counts();

      map_type x(ilist, 0, hasher(1), allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 2 * ilist.size());
      BOOST_TEST_EQ(raii::move_constructor, 2 * 11u);
    }
    check_raii_counts();
  }

  UNORDERED_AUTO_TEST (bucket_count_and_allocator) {
    raii::reset_counts();

    {
      map_type x(0, allocator_type(3));
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));
    }

    {
      map_type x(4096, allocator_type(3));
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));
    }
  }

  UNORDERED_AUTO_TEST (bucket_count_with_hasher_and_allocator) {
    raii::reset_counts();

    {
      map_type x(0, hasher(1), allocator_type(3));
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));
    }
  }

  template <class G>
  void iterator_range_with_bucket_count_and_allocator(
    G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    {
      allocator_type a(3);
      map_type x(values.begin(), values.end(), 0, a);
      test_fuzzy_matches_reference(x, reference_map, rg);

      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == a);
    }

    check_raii_counts();
  }

  template <class G>
  void iterator_range_with_bucket_count_hasher_and_allocator(
    G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    {
      allocator_type a(3);
      hasher hf(1);
      map_type x(values.begin(), values.end(), 0, hf, a);
      test_fuzzy_matches_reference(x, reference_map, rg);

      BOOST_TEST_EQ(x.hash_function(), hf);
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == a);
    }

    check_raii_counts();
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  from_iterator_range,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  copy_constructor,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  copy_constructor_with_insertion,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  move_constructor,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  move_constructor_with_insertion,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  iterator_range_with_allocator,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  iterator_range_with_bucket_count_and_allocator,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  iterator_range_with_bucket_count_hasher_and_allocator,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
