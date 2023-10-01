// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "exception_helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

using allocator_type = stateful_allocator<std::pair<raii const, raii> >;

using hasher = stateful_hash;
using key_equal = stateful_key_equal;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, allocator_type>;

namespace {
  test::seed_t initialize_seed(1794114520);

  template <class G> void copy_assign(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      unsigned num_throws = 0;

      auto begin = values.begin();
      auto mid =
        values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2);
      auto end = values.end();

      auto reference_map = boost::unordered_flat_map<raii, raii>(begin, mid);

      map_type x(
        begin, mid, values.size(), hasher(1), key_equal(2), allocator_type(3));

      map_type y(
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
      test_fuzzy_matches_reference(y, reference_map, rg);
    }
    check_raii_counts();
  }

  template <class G> void move_assign(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      unsigned num_throws = 0;

      auto begin = values.begin();
      auto mid =
        values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2);
      auto end = values.end();

      auto reference_map = boost::unordered_flat_map<raii, raii>(begin, mid);

      BOOST_TEST(
        !boost::allocator_is_always_equal<allocator_type>::type::value);

      BOOST_TEST(!boost::allocator_propagate_on_container_move_assignment<
                 allocator_type>::type::value);

      for (std::size_t i = 0; i < 2 * alloc_throw_threshold; ++i) {
        disable_exceptions();

        map_type x(begin, mid, values.size(), hasher(1), key_equal(2),
          allocator_type(3));

        map_type y(
          mid, end, values.size(), hasher(2), key_equal(1), allocator_type(4));

        enable_exceptions();
        try {
          y = std::move(x);
        } catch (...) {
          ++num_throws;
        }
        disable_exceptions();
        test_fuzzy_matches_reference(y, reference_map, rg);
      }

      BOOST_TEST_GT(num_throws, 0u);
    }
    check_raii_counts();
  }

  UNORDERED_AUTO_TEST (intializer_list_assign) {
    using value_type = typename map_type::value_type;

    std::initializer_list<value_type> values{
      value_type{raii{0}, raii{0}},
      value_type{raii{1}, raii{1}},
      value_type{raii{2}, raii{2}},
      value_type{raii{3}, raii{3}},
      value_type{raii{4}, raii{4}},
      value_type{raii{5}, raii{5}},
      value_type{raii{6}, raii{6}},
      value_type{raii{6}, raii{6}},
      value_type{raii{7}, raii{7}},
      value_type{raii{8}, raii{8}},
      value_type{raii{9}, raii{9}},
      value_type{raii{10}, raii{10}},
      value_type{raii{9}, raii{9}},
      value_type{raii{8}, raii{8}},
      value_type{raii{7}, raii{7}},
      value_type{raii{6}, raii{6}},
      value_type{raii{5}, raii{5}},
      value_type{raii{4}, raii{4}},
      value_type{raii{3}, raii{3}},
      value_type{raii{2}, raii{2}},
      value_type{raii{1}, raii{1}},
      value_type{raii{0}, raii{0}},
    };

    {
      raii::reset_counts();
      unsigned num_throws = 0;

      for (std::size_t i = 0; i < throw_threshold; ++i) {
        map_type x(0, hasher(1), key_equal(2), allocator_type(3));
        enable_exceptions();
        try {
          x = values;
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
  ((exception_value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  move_assign,
  ((exception_value_type_generator))
  ((default_generator)(sequential)))
// clang-format on

RUN_TESTS()
