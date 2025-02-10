// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "exception_helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

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

std::initializer_list<map_type::value_type> map_init_list{
  {raii{0}, raii{0}},
  {raii{1}, raii{1}},
  {raii{2}, raii{2}},
  {raii{3}, raii{3}},
  {raii{4}, raii{4}},
  {raii{5}, raii{5}},
  {raii{6}, raii{6}},
  {raii{6}, raii{6}},
  {raii{7}, raii{7}},
  {raii{8}, raii{8}},
  {raii{9}, raii{9}},
  {raii{10}, raii{10}},
  {raii{9}, raii{9}},
  {raii{8}, raii{8}},
  {raii{7}, raii{7}},
  {raii{6}, raii{6}},
  {raii{5}, raii{5}},
  {raii{4}, raii{4}},
  {raii{3}, raii{3}},
  {raii{2}, raii{2}},
  {raii{1}, raii{1}},
  {raii{0}, raii{0}},
};

std::initializer_list<set_type::value_type> set_init_list{
  raii{0},
  raii{1},
  raii{2},
  raii{3},
  raii{4},
  raii{5},
  raii{6},
  raii{6},
  raii{7},
  raii{8},
  raii{9},
  raii{10},
  raii{9},
  raii{8},
  raii{7},
  raii{6},
  raii{5},
  raii{4},
  raii{3},
  raii{2},
  raii{1},
  raii{0},
};

auto test_map_and_init_list=std::make_pair(test_map,map_init_list);
auto test_node_map_and_init_list=std::make_pair(test_node_map,map_init_list);
auto test_set_and_init_list=std::make_pair(test_set,set_init_list);
auto test_node_set_and_init_list=std::make_pair(test_node_set,set_init_list);

namespace {
  test::seed_t initialize_seed(795610904);

  template <class X>
  void bucket_constructor(X*)
  {
    raii::reset_counts();

    bool was_thrown = false;

    enable_exceptions();
    for (std::size_t i = 0; i < alloc_throw_threshold; ++i) {
      try {
        X m(128);
      } catch (...) {
        was_thrown = true;
      }
    }
    disable_exceptions();

    BOOST_TEST(was_thrown);
  }

  template <class X, class GF>
  void iterator_range(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      bool was_thrown = false;

      enable_exceptions();
      try {
        X x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
          allocator_type(3));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }

    {
      raii::reset_counts();

      bool was_thrown = false;

      enable_exceptions();
      try {
        X x(values.begin(), values.end(), allocator_type(3));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }

    {
      raii::reset_counts();

      bool was_thrown = false;

      enable_exceptions();
      try {
        X x(
          values.begin(), values.end(), values.size(), allocator_type(3));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }

    {
      raii::reset_counts();

      bool was_thrown = false;

      enable_exceptions();
      try {
        X x(values.begin(), values.end(), values.size(), hasher(1),
          allocator_type(3));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }
  }

  template <class X, class GF>
  void copy_constructor(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      bool was_thrown = false;

      try {
        X x(values.begin(), values.end(), 0);

        enable_exceptions();
        X y(x);
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }

    {
      raii::reset_counts();

      bool was_thrown = false;

      try {
        X x(values.begin(), values.end(), 0);

        enable_exceptions();
        X y(x, allocator_type(4));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }
  }

  template <class X, class GF>
  void move_constructor(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
 
    {
      raii::reset_counts();

      bool was_thrown = false;

      try {
        X x(values.begin(), values.end(), 0);

        enable_exceptions();
        X y(std::move(x), allocator_type(4));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }
  }

  template <class X, class IL>
  void initializer_list_bucket_count(std::pair<X*, IL> p)
  {
    using allocator_type = typename X::allocator_type;

    auto init_list = p.second;

    {
      raii::reset_counts();
      unsigned num_throws = 0;

      enable_exceptions();
      for (std::size_t i = 0; i < throw_threshold; ++i) {
        try {
          X x(init_list, 0, hasher(1), key_equal(2), allocator_type(3));
        } catch (...) {
          ++num_throws;
        }
      }
      disable_exceptions();

      BOOST_TEST_GT(num_throws, 0u);
      check_raii_counts();
    }

    {
      raii::reset_counts();
      unsigned num_throws = 0;

      enable_exceptions();
      for (std::size_t i = 0; i < alloc_throw_threshold * 2; ++i) {
        try {
          X x(init_list, allocator_type(3));
        } catch (...) {
          ++num_throws;
        }
      }
      disable_exceptions();

      BOOST_TEST_GT(num_throws, 0u);
      check_raii_counts();
    }

    {
      raii::reset_counts();
      unsigned num_throws = 0;

      enable_exceptions();
      for (std::size_t i = 0; i < alloc_throw_threshold * 2; ++i) {
        try {
          X x(init_list, init_list.size() * 2, allocator_type(3));
        } catch (...) {
          ++num_throws;
        }
      }
      disable_exceptions();

      BOOST_TEST_GT(num_throws, 0u);
      check_raii_counts();
    }

    {
      raii::reset_counts();
      unsigned num_throws = 0;

      enable_exceptions();
      for (std::size_t i = 0; i < throw_threshold; ++i) {
        try {
          X x(init_list, init_list.size() * 2, hasher(1), allocator_type(3));
        } catch (...) {
          ++num_throws;
        }
      }
      disable_exceptions();

      BOOST_TEST_GT(num_throws, 0u);
      check_raii_counts();
    }
  }
} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off
UNORDERED_TEST(
  bucket_constructor,
  ((test_map)(test_node_map)(test_set)(test_node_set)))

UNORDERED_TEST(
  iterator_range,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((exception_value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  copy_constructor,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((exception_value_type_generator_factory))
  ((default_generator)(sequential)))

UNORDERED_TEST(
  move_constructor,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((exception_value_type_generator_factory))
  ((default_generator)(sequential)))

UNORDERED_TEST(
  initializer_list_bucket_count,
  ((test_map_and_init_list)(test_node_map_and_init_list)
   (test_set_and_init_list)(test_node_set_and_init_list)))
// clang-format on

RUN_TESTS()
