// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "exception_helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

#include <boost/core/ignore_unused.hpp>

using allocator_type = stateful_allocator<std::pair<raii const, raii> >;

using hasher = stateful_hash;
using key_equal = stateful_key_equal;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, allocator_type>;

namespace {
  test::seed_t initialize_seed(223333016);

  template <class G> void merge(G gen, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    auto begin = values.begin();
    auto mid = begin + static_cast<std::ptrdiff_t>(values.size() / 2);
    auto end = values.end();

    {
      unsigned num_throws = 0;

      for (unsigned i = 0; i < 5 * alloc_throw_threshold; ++i) {
        disable_exceptions();

        map_type x1(0, hasher(1), key_equal(2), allocator_type(3));
        x1.insert(begin, mid);

        map_type x2(0, hasher(2), key_equal(1), allocator_type(3));
        x2.insert(mid, end);

        enable_exceptions();
        try {
          x1.merge(x2);
        } catch (...) {
          ++num_throws;
        }

        disable_exceptions();
        test_fuzzy_matches_reference(x1, reference_map, rg);
        test_fuzzy_matches_reference(x2, reference_map, rg);
      }

      BOOST_TEST_GT(num_throws, 0u);
    }

    check_raii_counts();
  }

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off
UNORDERED_TEST(
  merge,
  ((exception_value_type_generator))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
