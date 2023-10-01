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
  test::seed_t initialize_seed(795610904);

  UNORDERED_AUTO_TEST (bucket_constructor) {
    raii::reset_counts();

    bool was_thrown = false;

    enable_exceptions();
    for (std::size_t i = 0; i < alloc_throw_threshold; ++i) {
      try {
        map_type m(128);
      } catch (...) {
        was_thrown = true;
      }
    }
    disable_exceptions();

    BOOST_TEST(was_thrown);
  }

  template <class G> void iterator_range(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      bool was_thrown = false;

      enable_exceptions();
      try {
        map_type x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
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
        map_type x(values.begin(), values.end(), allocator_type(3));
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
        map_type x(
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
        map_type x(values.begin(), values.end(), values.size(), hasher(1),
          allocator_type(3));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }
  }

  template <class G> void copy_constructor(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      bool was_thrown = false;

      try {
        map_type x(values.begin(), values.end(), 0);

        enable_exceptions();
        map_type y(x);
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
        map_type x(values.begin(), values.end(), 0);

        enable_exceptions();
        map_type y(x, allocator_type(4));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }
  }

  template <class G> void move_constructor(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      bool was_thrown = false;

      try {
        map_type x(values.begin(), values.end(), 0);

        enable_exceptions();
        map_type y(std::move(x), allocator_type(4));
      } catch (...) {
        was_thrown = true;
      }
      disable_exceptions();

      BOOST_TEST(was_thrown);
      check_raii_counts();
    }
  }

  UNORDERED_AUTO_TEST (initializer_list_bucket_count) {
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

      enable_exceptions();
      for (std::size_t i = 0; i < throw_threshold; ++i) {
        try {
          map_type x(values, 0, hasher(1), key_equal(2), allocator_type(3));
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
          map_type x(values, allocator_type(3));
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
          map_type x(values, values.size() * 2, allocator_type(3));
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
          map_type x(values, values.size() * 2, hasher(1), allocator_type(3));
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
  iterator_range,
  ((exception_value_type_generator))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  copy_constructor,
  ((exception_value_type_generator))
  ((default_generator)(sequential)))

UNORDERED_TEST(
  move_constructor,
  ((exception_value_type_generator))
  ((default_generator)(sequential)))
// clang-format on

RUN_TESTS()
