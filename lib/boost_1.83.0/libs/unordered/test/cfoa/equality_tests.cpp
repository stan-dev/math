// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

test::seed_t initialize_seed{1634048962};

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

  UNORDERED_AUTO_TEST (simple_equality) {
    {
      map_type x1(
        {{1, 11}, {2, 22}}, 0, hasher(1), key_equal(2), allocator_type(3));

      map_type x2(
        {{1, 11}, {2, 22}}, 0, hasher(2), key_equal(2), allocator_type(3));

      map_type x3(
        {{1, 11}, {2, 23}}, 0, hasher(2), key_equal(2), allocator_type(3));

      map_type x4({{1, 11}}, 0, hasher(2), key_equal(2), allocator_type(3));

      BOOST_TEST_EQ(x1.size(), x2.size());
      BOOST_TEST(x1 == x2);
      BOOST_TEST(!(x1 != x2));

      BOOST_TEST_EQ(x1.size(), x3.size());
      BOOST_TEST(!(x1 == x3));
      BOOST_TEST(x1 != x3);

      BOOST_TEST(x1.size() != x4.size());
      BOOST_TEST(!(x1 == x4));
      BOOST_TEST(x1 != x4);
    }
  }

  template <class G> void insert_and_compare(G gen, test::random_generator rg)
  {
    auto vals1 = make_random_values(1024 * 8, [&] { return gen(rg); });
    boost::unordered_flat_map<raii, raii> reference_map(
      vals1.begin(), vals1.end());

    {
      raii::reset_counts();

      map_type x1(vals1.size(), hasher(1), key_equal(2), allocator_type(3));
      map_type x2(vals1.begin(), vals1.end(), vals1.size(), hasher(2),
        key_equal(2), allocator_type(3));

      std::thread t1, t2;

      std::mutex m;
      std::condition_variable cv;
      std::atomic_bool done{false};
      std::atomic<unsigned> num_compares{0};
      bool ready = false;

      BOOST_TEST(x1.empty());

      t1 = std::thread([&x1, &m, &cv, &vals1, &done, &ready] {
        for (std::size_t idx = 0; idx < vals1.size(); ++idx) {
          auto const& v = vals1[idx];
          x1.insert(v);

          if (idx % (vals1.size() / 128) == 0) {
            {
              std::unique_lock<std::mutex> lk(m);
              ready = true;
            }
            cv.notify_all();
          }
          std::this_thread::yield();
        }

        done = true;
        {
          std::unique_lock<std::mutex> lk(m);
          ready = true;
        }
        cv.notify_all();
      });

      t2 = std::thread([&x1, &x2, &m, &cv, &done, &num_compares, &ready] {
        do {
          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [&ready] { return ready; });
            ready = false;
          }

          volatile bool b = false;

          b = x1 == x2;
          b = x1 != x2;

          b;

          ++num_compares;
          std::this_thread::yield();
        } while (!done);

        BOOST_TEST(done);
      });

      t1.join();
      t2.join();

      BOOST_TEST_GE(num_compares, 1u);

      BOOST_TEST(x1 == x2);
      BOOST_TEST(!(x1 != x2));

      test_matches_reference(x1, reference_map);
    }
    check_raii_counts();
  }
} // namespace

// clang-format off
UNORDERED_TEST(
  insert_and_compare,
  ((value_type_generator))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
