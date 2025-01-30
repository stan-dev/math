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
  test::seed_t initialize_seed(1794114520);

  template <class X, class GF>
  void copy_assign(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      unsigned num_throws = 0;

      auto begin = values.begin();
      auto mid =
        values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2);
      auto end = values.end();

      auto reference_cont = reference_container<X>(begin, mid);

      X x(
        begin, mid, values.size(), hasher(1), key_equal(2), allocator_type(3));

      X y(
        mid, end, values.size(), hasher(2), key_equal(1), allocator_type(4));

      BOOST_TEST(!y.empty());

      enable_exceptions();
      for (std::size_t i = 0; i < 2 * alloc_throw_threshold; ++i) {
        try {
          y = x;
        } catch (...) {
          ++num_throws;
        }
      }

      disable_exceptions();

      BOOST_TEST_GT(num_throws, 0u);
      test_fuzzy_matches_reference(y, reference_cont, rg);
    }
    check_raii_counts();
  }

  template <class X, class GF>
  void move_assign(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      unsigned num_throws = 0;

      auto begin = values.begin();
      auto mid =
        values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2);
      auto end = values.end();

      auto reference_cont = reference_container<X>(begin, mid);

      BOOST_TEST(
        !boost::allocator_is_always_equal<allocator_type>::type::value);

      BOOST_TEST(!boost::allocator_propagate_on_container_move_assignment<
                 allocator_type>::type::value);

      for (std::size_t i = 0; i < 2 * alloc_throw_threshold; ++i) {
        disable_exceptions();

        X x(begin, mid, values.size(), hasher(1), key_equal(2),
          allocator_type(3));

        X y(
          mid, end, values.size(), hasher(2), key_equal(1), allocator_type(4));

        enable_exceptions();
        try {
          y = std::move(x);
        } catch (...) {
          ++num_throws;
        }
        disable_exceptions();
        test_fuzzy_matches_reference(y, reference_cont, rg);
      }

      BOOST_TEST_GT(num_throws, 0u);
    }
    check_raii_counts();
  }

  template <class X, class IL>
  void intializer_list_assign(std::pair<X*, IL> p)
  {
    using allocator_type = typename X::allocator_type;

    auto init_list = p.second;

    {
      raii::reset_counts();
      unsigned num_throws = 0;

      for (std::size_t i = 0; i < throw_threshold; ++i) {
        X x(0, hasher(1), key_equal(2), allocator_type(3));
        enable_exceptions();
        try {
          x = init_list;
        } catch (...) {
          ++num_throws;
        }
        disable_exceptions();
      }

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
  copy_assign,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((exception_value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  move_assign,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((exception_value_type_generator_factory))
  ((default_generator)(sequential)))

UNORDERED_TEST(
  intializer_list_assign,
  ((test_map_and_init_list)(test_node_map_and_init_list)
   (test_set_and_init_list)(test_node_set_and_init_list)))
// clang-format on

RUN_TESTS()
