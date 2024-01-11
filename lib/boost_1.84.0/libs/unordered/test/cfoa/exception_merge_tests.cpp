// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "exception_helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>

#include <boost/core/ignore_unused.hpp>

using hasher = stateful_hash;
using key_equal = stateful_key_equal;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, stateful_allocator<std::pair<raii const, raii> > >;

using set_type = boost::unordered::concurrent_flat_set<raii, hasher,
  key_equal, stateful_allocator<raii> >;

map_type* test_map;
set_type* test_set;

namespace {
  test::seed_t initialize_seed(223333016);

  template <class X, class GF>
  void merge(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    auto begin = values.begin();
    auto mid = begin + static_cast<std::ptrdiff_t>(values.size() / 2);
    auto end = values.end();

    {
      unsigned num_throws = 0;

      for (unsigned i = 0; i < 5 * alloc_throw_threshold; ++i) {
        disable_exceptions();

        X x1(0, hasher(1), key_equal(2), allocator_type(3));
        x1.insert(begin, mid);

        X x2(0, hasher(2), key_equal(1), allocator_type(3));
        x2.insert(mid, end);

        enable_exceptions();
        try {
          x1.merge(x2);
        } catch (...) {
          ++num_throws;
        }

        disable_exceptions();
        test_fuzzy_matches_reference(x1, reference_cont, rg);
        test_fuzzy_matches_reference(x2, reference_cont, rg);
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
  ((test_map)(test_set))
  ((exception_value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
