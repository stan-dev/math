// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

test::seed_t initialize_seed{402031699};

using test::default_generator;
using test::limited_range;
using test::sequential;

using hasher = stateful_hash;
using key_equal = stateful_key_equal;
using allocator_type = stateful_allocator<std::pair<raii const, raii> >;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, allocator_type>;

using map_value_type = typename map_type::value_type;

struct
{
  template <class X1, class X2>
  std::size_t operator()(X1& x1, X2& x2) const noexcept
  {
    return x1.merge(x2);
  }
} lvalue_merge;

struct
{
  template <class X1, class X2>
  std::size_t operator()(X1& x1, X2& x2) const noexcept
  {
    return x1.merge(std::move(x2));
  }
} rvalue_merge;

namespace {
  template <class F, class G>
  void merge_tests(F merger, G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 8, [&] { return gen(rg); });

    auto ref_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    {
      raii::reset_counts();

      map_type x(values.size(), hasher(1), key_equal(2), allocator_type(3));

      auto const old_cc = +raii::copy_constructor;

      std::atomic<unsigned long long> expected_copies{0};
      std::atomic<unsigned long long> num_merged{0};

      thread_runner(values, [&x, &expected_copies, &num_merged, merger](
                              boost::span<map_value_type> s) {
        using map2_type = boost::unordered::concurrent_flat_map<raii, raii,
          std::hash<raii>, std::equal_to<raii>, allocator_type>;

        map2_type y(s.size(), allocator_type(3));
        for (auto const& v : s) {
          y.insert(v);
        }
        expected_copies += 2 * y.size();

        BOOST_TEST(x.get_allocator() == y.get_allocator());
        num_merged += merger(x, y);
      });

      BOOST_TEST_EQ(raii::copy_constructor, old_cc + expected_copies);
      BOOST_TEST_EQ(raii::move_constructor, 2 * ref_map.size());
      BOOST_TEST_EQ(+num_merged, ref_map.size());

      test_fuzzy_matches_reference(x, ref_map, rg);
    }
    check_raii_counts();
  }

  template <class G>
  void insert_and_merge_tests(G gen, test::random_generator rg)
  {
    using map2_type = boost::unordered::concurrent_flat_map<raii, raii,
      std::hash<raii>, std::equal_to<raii>, allocator_type>;

    auto vals1 = make_random_values(1024 * 8, [&] { return gen(rg); });
    auto vals2 = make_random_values(1024 * 4, [&] { return gen(rg); });

    auto ref_map = boost::unordered_flat_map<raii, raii>();
    ref_map.insert(vals1.begin(), vals1.end());
    ref_map.insert(vals2.begin(), vals2.end());

    {
      raii::reset_counts();

      map_type x1(2 * vals1.size(), hasher(1), key_equal(2), allocator_type(3));

      map2_type x2(2 * vals1.size(), allocator_type(3));

      std::thread t1, t2, t3;
      boost::compat::latch l(2);

      std::mutex m;
      std::condition_variable cv;
      std::atomic_bool done1{false}, done2{false};
      std::atomic<unsigned long long> num_merges{0};
      std::atomic<unsigned long long> call_count{0};
      bool ready = false;

      auto const old_mc = +raii::move_constructor;
      BOOST_TEST_EQ(old_mc, 0u);

      t1 = std::thread([&x1, &vals1, &l, &done1, &cv, &ready, &m] {
        l.arrive_and_wait();

        for (std::size_t idx = 0; idx < vals1.size(); ++idx) {
          auto const& val = vals1[idx];
          x1.insert(val);

          if (idx % (vals1.size() / 128) == 0) {
            {
              std::unique_lock<std::mutex> lk(m);
              ready = true;
            }
            cv.notify_all();
            std::this_thread::yield();
          }
        }

        done1 = true;
        {
          std::unique_lock<std::mutex> lk(m);
          ready = true;
        }
        cv.notify_all();
      });

      t2 = std::thread([&x2, &vals2, &l, &done2, &cv, &m, &ready] {
        l.arrive_and_wait();

        for (std::size_t idx = 0; idx < vals2.size(); ++idx) {
          auto const& val = vals2[idx];
          x2.insert(val);
          if (idx % 100 == 0) {
            std::this_thread::yield();
          }
        }

        done2 = true;
        {
          std::unique_lock<std::mutex> lk(m);
          ready = true;
        }
        cv.notify_all();
      });

      t3 = std::thread(
        [&x1, &x2, &m, &cv, &done1, &done2, &num_merges, &call_count, &ready] {
          while (x1.empty() && x2.empty()) {
          }

          do {
            {
              std::unique_lock<std::mutex> lk(m);
              cv.wait(lk, [&ready] { return ready; });
              ready = false;
            }

            num_merges += x1.merge(x2);
            std::this_thread::yield();
            num_merges += x2.merge(x1);

            call_count += 1;

          } while (!done1 || !done2);

          BOOST_TEST(done1);
          BOOST_TEST(done2);
        });

      t1.join();
      t2.join();
      t3.join();

      if (num_merges > 0) {
        // num merges is 0 most commonly in the cast of the limited_range
        // generator as both maps will contains keys from 0 to 99
        BOOST_TEST_EQ(+raii::move_constructor, 2 * num_merges);
        BOOST_TEST_GE(call_count, 1u);
      }

      x1.merge(x2);
      test_fuzzy_matches_reference(x1, ref_map, rg);
    }

    check_raii_counts();
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  merge_tests,
  ((lvalue_merge)(rvalue_merge))
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  insert_and_merge_tests,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
