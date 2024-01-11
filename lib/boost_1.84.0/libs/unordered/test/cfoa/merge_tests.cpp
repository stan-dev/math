// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>

test::seed_t initialize_seed{402031699};

using test::default_generator;
using test::limited_range;
using test::sequential;

using hasher = stateful_hash;
using key_equal = stateful_key_equal;

using map_type = boost::unordered::concurrent_flat_map<raii, raii,
  hasher, key_equal, stateful_allocator<std::pair<raii const, raii> > >;
using map2_type = boost::unordered::concurrent_flat_map<raii, raii,
  std::hash<raii>, std::equal_to<raii>,
  stateful_allocator<std::pair<raii const, raii> > >;

using set_type = boost::unordered::concurrent_flat_set<raii, hasher,
  key_equal, stateful_allocator<raii> >;
using set2_type = boost::unordered::concurrent_flat_set<raii, std::hash<raii>,
  std::equal_to<raii>, stateful_allocator<raii> >;

map_type* test_map;
map2_type* test_map2;
auto test_maps=std::make_pair(test_map,test_map2);

set_type* test_set;
set2_type* test_set2;
auto test_sets=std::make_pair(test_set,test_set2);

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
  template <typename X, typename Y, class F, class GF>
  void merge_tests(
    std::pair<X*, Y*>, F merger, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 8, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));

      auto const old_cc = +raii::copy_constructor;

      std::atomic<unsigned long long> expected_copies{0};
      std::atomic<unsigned long long> num_merged{0};

      thread_runner(values, [&x, &expected_copies, &num_merged, merger](
                              boost::span<value_type> s) {
        Y y(s.size(), allocator_type(3));
        for (auto const& v : s) {
          y.insert(v);
        }
        expected_copies += value_type_cardinality * y.size();

        BOOST_TEST(x.get_allocator() == y.get_allocator());
        num_merged += merger(x, y);
      });

      BOOST_TEST_EQ(raii::copy_constructor, old_cc + expected_copies);
      BOOST_TEST_EQ(
        raii::move_constructor, 
        value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(+num_merged, reference_cont.size());

      test_fuzzy_matches_reference(x, reference_cont, rg);
    }
    check_raii_counts();
  }

  template <typename X, typename Y, class GF>
  void insert_and_merge_tests(
    std::pair<X*, Y*>, GF gen_factory, test::random_generator rg)
  {
    static constexpr auto value_type_cardinality = 
      value_cardinality<typename X::value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto vals1 = make_random_values(1024 * 8, [&] { return gen(rg); });
    auto vals2 = make_random_values(1024 * 4, [&] { return gen(rg); });

    auto reference_cont = reference_container<X>();
    reference_cont.insert(vals1.begin(), vals1.end());
    reference_cont.insert(vals2.begin(), vals2.end());

    {
      raii::reset_counts();

      X x1(2 * vals1.size(), hasher(1), key_equal(2), allocator_type(3));

      Y x2(2 * vals1.size(), allocator_type(3));

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
        BOOST_TEST_EQ(
          +raii::move_constructor, value_type_cardinality * num_merges);
        BOOST_TEST_GE(call_count, 1u);
      }

      x1.merge(x2);
      test_fuzzy_matches_reference(x1, reference_cont, rg);
    }

    check_raii_counts();
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  merge_tests,
  ((test_maps)(test_sets))
  ((lvalue_merge)(rvalue_merge))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  insert_and_merge_tests,
  ((test_maps)(test_sets))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
