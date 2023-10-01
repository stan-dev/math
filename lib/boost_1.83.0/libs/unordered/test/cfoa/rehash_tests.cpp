// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

using test::default_generator;
using test::limited_range;
using test::sequential;

using hasher = stateful_hash;
using key_equal = stateful_key_equal;
using allocator_type = stateful_allocator<std::pair<raii const, raii> >;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, allocator_type>;

using map_value_type = typename map_type::value_type;

namespace {
  test::seed_t initialize_seed{748775921};

  UNORDERED_AUTO_TEST (rehash_no_insert) {
    map_type x(0, hasher(1), key_equal(2), allocator_type(3));
    BOOST_TEST_EQ(x.bucket_count(), 0u);

    x.rehash(1024);
    BOOST_TEST_GE(x.bucket_count(), 1024u);

    x.rehash(512);
    BOOST_TEST_GE(x.bucket_count(), 512u);
    BOOST_TEST_LT(x.bucket_count(), 1024u);

    x.rehash(0);
    BOOST_TEST_EQ(x.bucket_count(), 0u);
  }

  UNORDERED_AUTO_TEST (reserve_no_insert) {
    using size_type = map_type::size_type;

    map_type x(0, hasher(1), key_equal(2), allocator_type(3));

    auto f = [&x](double c) {
      return static_cast<size_type>(std::ceil(c / x.max_load_factor()));
    };

    BOOST_TEST_EQ(x.bucket_count(), f(0.0));

    x.reserve(1024);
    BOOST_TEST_GE(x.bucket_count(), f(1024.0));

    x.reserve(512);
    BOOST_TEST_GE(x.bucket_count(), f(512.0));
    BOOST_TEST_LT(x.bucket_count(), f(1024.0));

    x.reserve(0);
    BOOST_TEST_EQ(x.bucket_count(), f(0.0));
  }

  template <class G>
  void insert_and_erase_with_rehash(G gen, test::random_generator rg)
  {
    auto vals1 = make_random_values(1024 * 8, [&] { return gen(rg); });

    auto erase_indices = std::vector<std::size_t>(vals1.size());
    for (std::size_t idx = 0; idx < erase_indices.size(); ++idx) {
      erase_indices[idx] = idx;
    }
    shuffle_values(erase_indices);

    auto ref_map = boost::unordered_flat_map<raii, raii>();
    ref_map.insert(vals1.begin(), vals1.end());

    {
      raii::reset_counts();

      map_type x(0, hasher(1), key_equal(2), allocator_type(3));

      std::thread t1, t2, t3;
      boost::compat::latch l(2);

      std::mutex m;
      std::condition_variable cv;
      std::atomic_bool done1{false}, done2{false};
      std::atomic<unsigned long long> call_count{0};
      bool ready = false;

      auto const old_mc = +raii::move_constructor;
      BOOST_TEST_EQ(old_mc, 0u);

      t1 = std::thread([&x, &vals1, &l, &done1, &cv, &ready, &m] {
        l.arrive_and_wait();

        for (std::size_t idx = 0; idx < vals1.size(); ++idx) {
          auto const& val = vals1[idx];
          x.insert(val);

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

      t2 =
        std::thread([&x, &vals1, &erase_indices, &l, &done2, &cv, &m, &ready] {
          l.arrive_and_wait();

          for (std::size_t idx = 0; idx < erase_indices.size(); ++idx) {
            auto const& val = vals1[erase_indices[idx]];
            x.erase(val.first);
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

      t3 =
        std::thread([&x, &vals1, &m, &cv, &done1, &done2, &call_count, &ready] {
          do {
            {
              std::unique_lock<std::mutex> lk(m);
              cv.wait(lk, [&ready] { return ready; });
              ready = false;
            }

            auto const bc = static_cast<std::size_t>(rand()) % vals1.size();
            x.rehash(bc);
            call_count += 1;

            std::this_thread::yield();
          } while (!done1 || !done2);

          BOOST_TEST(done1);
          BOOST_TEST(done2);
        });

      t1.join();
      t2.join();
      t3.join();

      BOOST_TEST_GE(call_count, 1u);

      test_fuzzy_matches_reference(x, ref_map, rg);
    }

    check_raii_counts();
  }
} // namespace

// clang-format off
UNORDERED_TEST(
  insert_and_erase_with_rehash,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
