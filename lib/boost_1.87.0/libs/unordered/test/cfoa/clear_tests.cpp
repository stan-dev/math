// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

test::seed_t initialize_seed{674140082};

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
  template <class X, class GF>
  void clear_tests(X*, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    raii::reset_counts();

    X x(values.begin(), values.end(), values.size(), hasher(1),
      key_equal(2), allocator_type(3));

    auto const old_size = x.size();
    auto const old_d = +raii::destructor;

    thread_runner(values, [&x](boost::span<value_type> s) {
      (void)s;
      x.clear();
    });

    BOOST_TEST(x.empty());
    BOOST_TEST_EQ(raii::destructor, old_d + value_type_cardinality * old_size);

    check_raii_counts();
  }

  template <class X, class GF>
  void insert_and_clear(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    std::thread t1, t2;

    {
      X x(0, hasher(1), key_equal(2), allocator_type(3));

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
        test_fuzzy_matches_reference(x, reference_cont, rg);
      }
    }

    check_raii_counts();
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  clear_tests,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(insert_and_clear,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
