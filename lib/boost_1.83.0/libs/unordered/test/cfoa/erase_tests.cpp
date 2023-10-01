// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

#include <boost/core/ignore_unused.hpp>

namespace {
  test::seed_t initialize_seed(3292023);

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

      BOOST_TEST_EQ(raii::default_constructor + raii::copy_constructor +
                      raii::move_constructor,
        raii::destructor + 2 * x.size());

      thread_runner(values, [&values, &num_erased, &x](boost::span<T>) {
        for (auto const& k : values) {
          auto count = x.erase(k.first);
          num_erased += count;
          BOOST_TEST_LE(count, 1u);
          BOOST_TEST_GE(count, 0u);
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * old_size);

      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(num_erased, old_size);
    }
  } lvalue_eraser;

  struct transp_lvalue_eraser_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_erased{0};
      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      BOOST_TEST_EQ(raii::default_constructor + raii::copy_constructor +
                      raii::move_constructor,
        raii::destructor + 2 * x.size());

      thread_runner(values, [&num_erased, &x](boost::span<T> s) {
        for (auto const& k : s) {
          auto count = x.erase(k.first.x_);
          num_erased += count;
          BOOST_TEST_LE(count, 1u);
          BOOST_TEST_GE(count, 0u);
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * num_erased);

      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(num_erased, old_size);
    }
  } transp_lvalue_eraser;

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

      thread_runner(values, [&num_erased, &x, threshold](boost::span<T> s) {
        for (auto const& k : s) {
          auto count = x.erase_if(k.first,
            [threshold](value_type& v) { return v.second.x_ > threshold; });
          num_erased += count;
          BOOST_TEST_LE(count, 1u);
          BOOST_TEST_GE(count, 0u);
        }
      });

      BOOST_TEST_EQ(num_erased, expected_erasures);
      BOOST_TEST_EQ(x.size(), old_size - num_erased);

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * num_erased);
    }
  } lvalue_eraser_if;

  struct transp_lvalue_eraser_if_type
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

      thread_runner(values, [&num_erased, &x, threshold](boost::span<T> s) {
        for (auto const& k : s) {
          auto count = x.erase_if(k.first.x_,
            [threshold](value_type& v) { return v.second.x_ > threshold; });
          num_erased += count;
          BOOST_TEST_LE(count, 1u);
          BOOST_TEST_GE(count, 0u);
        }
      });

      BOOST_TEST_EQ(num_erased, expected_erasures);
      BOOST_TEST_EQ(x.size(), old_size - num_erased);

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * num_erased);
    }
  } transp_lvalue_eraser_if;

  struct erase_if_type
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

      thread_runner(
        values, [&num_erased, &x, threshold](boost::span<T> /* s */) {
          for (std::size_t i = 0; i < 128; ++i) {
            auto count = x.erase_if(
              [threshold](value_type& v) { return v.second.x_ > threshold; });
            num_erased += count;
          }
        });

      BOOST_TEST_EQ(num_erased, expected_erasures);
      BOOST_TEST_EQ(x.size(), old_size - num_erased);

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * num_erased);
    }
  } erase_if;

  struct free_fn_erase_if_type
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

      thread_runner(
        values, [&num_erased, &x, threshold](boost::span<T> /* s */) {
          for (std::size_t i = 0; i < 128; ++i) {
            auto count = boost::unordered::erase_if(x,
              [threshold](value_type& v) { return v.second.x_ > threshold; });
            num_erased += count;
          }
        });

      BOOST_TEST_EQ(num_erased, expected_erasures);
      BOOST_TEST_EQ(x.size(), old_size - num_erased);

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * num_erased);
    }
  } free_fn_erase_if;

  struct erase_if_exec_policy_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
#if defined(BOOST_UNORDERED_PARALLEL_ALGORITHMS)
      using value_type = typename X::value_type;

      std::atomic<std::uint64_t> num_invokes{0};

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

      thread_runner(values, [&num_invokes, &x, threshold](boost::span<T> s) {
        (void)s;
        x.erase_if(
          std::execution::par, [&num_invokes, threshold](value_type& v) {
            ++num_invokes;
            return v.second.x_ > threshold;
          });
      });

      BOOST_TEST_GE(+num_invokes, old_size);
      BOOST_TEST_LE(+num_invokes, old_size * num_threads);

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(raii::destructor, old_d + 2 * expected_erasures);
#else
      (void)values;
      (void)x;
#endif
    }
  } erase_if_exec_policy;

  template <class X, class G, class F>
  void erase(X*, G gen, F eraser, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());
    raii::reset_counts();

    {
      X x;

      x.insert(values.begin(), values.end());

      BOOST_TEST_EQ(x.size(), reference_map.size());

      test_fuzzy_matches_reference(x, reference_map, rg);

      eraser(values, x);
      test_fuzzy_matches_reference(x, reference_map, rg);
    }

    check_raii_counts();
  }

  boost::unordered::concurrent_flat_map<raii, raii>* map;
  boost::unordered::concurrent_flat_map<raii, raii, transp_hash,
    transp_key_equal>* transparent_map;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off
UNORDERED_TEST(
  erase,
  ((map))
  ((value_type_generator)(init_type_generator))
  ((lvalue_eraser)(lvalue_eraser_if)(erase_if)(free_fn_erase_if)(erase_if_exec_policy))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  erase,
  ((transparent_map))
  ((value_type_generator)(init_type_generator))
  ((transp_lvalue_eraser)(transp_lvalue_eraser_if)(erase_if_exec_policy))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
