// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "exception_helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

#include <boost/core/ignore_unused.hpp>

namespace {
  test::seed_t initialize_seed(3202923);

  struct lvalue_eraser_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_erased{0};

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      enable_exceptions();
      thread_runner(values, [&values, &num_erased, &x](boost::span<T>) {
        for (auto const& k : values) {
          try {
            auto count = x.erase(k.first);
            BOOST_TEST_LE(count, 1u);
            BOOST_TEST_GE(count, 0u);

            num_erased += count;
          } catch (...) {
          }
        }
      });
      disable_exceptions();

      BOOST_TEST_EQ(x.size(), old_size - num_erased);

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * num_erased);
    }
  } lvalue_eraser;

  struct lvalue_eraser_if_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using value_type = typename X::value_type;

      std::atomic<std::uint64_t> num_erased{0};

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      auto max = 0;
      x.visit_all([&max](value_type const& v) {
        if (v.second.x_ > max) {
          max = v.second.x_;
        }
      });

      auto threshold = max / 2;

      auto expected_erasures = 0u;
      x.visit_all([&expected_erasures, threshold](value_type const& v) {
        if (v.second.x_ > threshold) {
          ++expected_erasures;
        }
      });

      enable_exceptions();
      thread_runner(values, [&num_erased, &x, threshold](boost::span<T> s) {
        for (auto const& k : s) {
          try {
            auto count = x.erase_if(k.first,
              [threshold](value_type& v) { return v.second.x_ > threshold; });
            num_erased += count;
            BOOST_TEST_LE(count, 1u);
            BOOST_TEST_GE(count, 0u);
          } catch (...) {
          }
        }
      });
      disable_exceptions();

      BOOST_TEST_LE(num_erased, expected_erasures);
      BOOST_TEST_EQ(x.size(), old_size - num_erased);

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * num_erased);
    }
  } lvalue_eraser_if;

  struct erase_if_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using value_type = typename X::value_type;

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      auto max = 0;
      x.visit_all([&max](value_type const& v) {
        if (v.second.x_ > max) {
          max = v.second.x_;
        }
      });

      auto threshold = max / 2;

      auto expected_erasures = 0u;
      x.visit_all([&expected_erasures, threshold](value_type const& v) {
        if (v.second.x_ > threshold) {
          ++expected_erasures;
        }
      });

      enable_exceptions();
      thread_runner(values, [&x, threshold](boost::span<T> /* s */) {
        for (std::size_t i = 0; i < 256; ++i) {
          try {
            x.erase_if([threshold](value_type& v) {
              static std::atomic<std::uint32_t> c{0};
              auto t = ++c;
              if (should_throw && (t % throw_threshold == 0)) {
                throw exception_tag{};
              }

              return v.second.x_ > threshold;
            });
          } catch (...) {
          }
        }
      });
      disable_exceptions();

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * (old_size - x.size()));
    }
  } erase_if;

  struct free_fn_erase_if_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using value_type = typename X::value_type;

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      auto max = 0;
      x.visit_all([&max](value_type const& v) {
        if (v.second.x_ > max) {
          max = v.second.x_;
        }
      });

      auto threshold = max / 2;

      enable_exceptions();
      thread_runner(values, [&x, threshold](boost::span<T> /* s */) {
        for (std::size_t i = 0; i < 256; ++i) {
          try {
            boost::unordered::erase_if(x, [threshold](value_type& v) {
              static std::atomic<std::uint32_t> c{0};
              auto t = ++c;
              if (should_throw && (t % throw_threshold == 0)) {
                throw exception_tag{};
              }

              return v.second.x_ > threshold;
            });

          } catch (...) {
          }
        }
      });
      disable_exceptions();

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * (old_size - x.size()));
    }
  } free_fn_erase_if;

  template <class X, class G, class F>
  void erase(X*, G gen, F eraser, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    {
      X x(values.size());
      for (auto const& v : values) {
        x.insert(v);
      }

      BOOST_TEST_EQ(x.size(), reference_map.size());
      BOOST_TEST_EQ(raii::destructor, 0u);

      test_fuzzy_matches_reference(x, reference_map, rg);

      eraser(values, x);
      test_fuzzy_matches_reference(x, reference_map, rg);
    }

    check_raii_counts();
  }

  boost::unordered::concurrent_flat_map<raii, raii, stateful_hash,
    stateful_key_equal, stateful_allocator<std::pair<raii const, raii> > >* map;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off
UNORDERED_TEST(
  erase,
  ((map))
  ((exception_value_type_generator)(exception_init_type_generator))
  ((lvalue_eraser)(lvalue_eraser_if)(erase_if)(free_fn_erase_if))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
