// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

test::seed_t initialize_seed{674140082};

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
  template <class G> void clear_tests(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    raii::reset_counts();

    map_type x(values.begin(), values.end(), values.size(), hasher(1),
      key_equal(2), allocator_type(3));

    auto const old_size = x.size();
    auto const old_d = +raii::destructor;

    thread_runner(values, [&x](boost::span<map_value_type> s) {
      (void)s;
      x.clear();
    });

    BOOST_TEST(x.empty());
    BOOST_TEST_EQ(raii::destructor, old_d + 2 * old_size);

    check_raii_counts();
  }

  template <class G> void insert_and_clear(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    std::thread t1, t2;

    {
      map_type x(0, hasher(1), key_equal(2), allocator_type(3));

      std::mutex m;
      std::condition_variable cv;
      std::atomic<bool> done{false};
      std::atomic<unsigned> num_clears{0};

      bool ready = false;

      t1 = std::thread([&x, &values, &cv, &done, &m, &ready] {
        for (auto i = 0u; i < values.size(); ++i) {
          x.insert(values[i]);
          if (i % (values.size() / 128) == 0) {
            {
              std::unique_lock<std::mutex> lk(m);
              ready = true;
            }
            cv.notify_all();
          }
        }

        done = true;
        {
          std::unique_lock<std::mutex> lk(m);
          ready = true;
        }
        cv.notify_all();
      });

      t2 = std::thread([&x, &m, &cv, &done, &ready, &num_clears] {
        do {
          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [&ready] { return ready; });
            ready = false;
          }
          x.clear();
          ++num_clears;
        } while (!done);
      });

      t1.join();
      t2.join();

      BOOST_TEST_GE(num_clears, 1u);

      if (!x.empty()) {
        test_fuzzy_matches_reference(x, reference_map, rg);
      }
    }

    check_raii_counts();
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  clear_tests,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(insert_and_clear,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
