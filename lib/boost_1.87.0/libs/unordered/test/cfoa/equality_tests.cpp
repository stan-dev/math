// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

test::seed_t initialize_seed{1634048962};

using test::default_generator;
using test::limited_range;
using test::sequential;

using hasher = stateful_hash;
using key_equal = stateful_key_equal;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, stateful_allocator<std::pair<raii const, raii> > >;

using node_map_type = boost::unordered::concurrent_node_map<raii, raii, hasher,
  key_equal, stateful_allocator<std::pair<raii const, raii> > >;

using set_type = boost::unordered::concurrent_flat_set<raii, hasher,
  key_equal, stateful_allocator<raii> >;

using node_set_type = boost::unordered::concurrent_node_set<raii, hasher,
  key_equal, stateful_allocator<raii> >;

map_type* test_map;
node_map_type* test_node_map;
set_type* test_set;
node_set_type* test_node_set;

namespace {

  template <class X>
  void simple_map_equality(X*)
  {
    using allocator_type = typename X::allocator_type;

    {
      X x1(
        {{1, 11}, {2, 22}}, 0, hasher(1), key_equal(2), allocator_type(3));

      X x2(
        {{1, 11}, {2, 22}}, 0, hasher(2), key_equal(2), allocator_type(3));

      X x3(
        {{1, 11}, {2, 23}}, 0, hasher(2), key_equal(2), allocator_type(3));

      X x4({{1, 11}}, 0, hasher(2), key_equal(2), allocator_type(3));

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

  template <class X>
  void simple_set_equality(X*)
  {
    using allocator_type = typename X::allocator_type;

    {
      X x1(
        {1, 2}, 0, hasher(1), key_equal(2), allocator_type(3));

      X x2(
        {1, 2}, 0, hasher(2), key_equal(2), allocator_type(3));

      X x3({1}, 0, hasher(2), key_equal(2), allocator_type(3));

      BOOST_TEST_EQ(x1.size(), x2.size());
      BOOST_TEST(x1 == x2);
      BOOST_TEST(!(x1 != x2));

      BOOST_TEST(x1.size() != x3.size());
      BOOST_TEST(!(x1 == x3));
      BOOST_TEST(x1 != x3);
    }
  }

  template <class X, class GF>
  void insert_and_compare(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto vals1 = make_random_values(1024 * 8, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(vals1.begin(), vals1.end());

    {
      raii::reset_counts();

      X x1(vals1.size(), hasher(1), key_equal(2), allocator_type(3));
      X x2(vals1.begin(), vals1.end(), vals1.size(), hasher(2),
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

      test_matches_reference(x1, reference_cont);
    }
    check_raii_counts();
  }
} // namespace

// clang-format off
UNORDERED_TEST(
  simple_map_equality,
  ((test_map)(test_node_map)))

UNORDERED_TEST(
  simple_set_equality,
  ((test_set)(test_node_set)))

UNORDERED_TEST(
  insert_and_compare,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
