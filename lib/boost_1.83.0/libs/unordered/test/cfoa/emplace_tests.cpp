// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

#include <boost/core/ignore_unused.hpp>

namespace {
  test::seed_t initialize_seed(335740237);

  struct lvalue_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto const& r : s) {
          bool b = x.emplace(r.first.x_, r.second.x_);
          if (b) {
            ++num_inserts;
          }
        }
      });
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(raii::default_constructor, 2 * values.size());

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, 2 * x.size());

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_emplacer;

  struct norehash_lvalue_emplacer_type : public lvalue_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      x.reserve(values.size());
      lvalue_emplacer_type::operator()(values, x);
      BOOST_TEST_EQ(raii::move_constructor, 2 * x.size());
    }
  } norehash_lvalue_emplacer;

  struct lvalue_emplace_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.emplace_or_cvisit(
            r.first.x_, r.second.x_,
            [&num_invokes](typename X::value_type const& v) {
              (void)v;
              ++num_invokes;
            });

          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 2 * values.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, 2 * x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_emplace_or_cvisit;

  struct lvalue_emplace_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.emplace_or_visit(
            r.first.x_, r.second.x_,
            [&num_invokes](typename X::value_type& v) {
              (void)v;
              ++num_invokes;
            });

          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 2 * values.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, 2 * x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_emplace_or_visit;

  template <class X, class G, class F>
  void emplace(X*, G gen, F emplacer, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());
    raii::reset_counts();

    {
      X x;

      emplacer(values, x);

      BOOST_TEST_EQ(x.size(), reference_map.size());

      using value_type = typename X::value_type;
      BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& kv) {
        BOOST_TEST(reference_map.contains(kv.first));
        if (rg == test::sequential) {
          BOOST_TEST_EQ(kv.second, reference_map[kv.first]);
        }
      }));
    }

    BOOST_TEST_GE(raii::default_constructor, 0u);
    BOOST_TEST_GE(raii::copy_constructor, 0u);
    BOOST_TEST_GE(raii::move_constructor, 0u);
    BOOST_TEST_GT(raii::destructor, 0u);

    BOOST_TEST_EQ(raii::default_constructor + raii::copy_constructor +
                    raii::move_constructor,
      raii::destructor);
  }

  boost::unordered::concurrent_flat_map<raii, raii>* map;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off

UNORDERED_TEST(
  emplace,
  ((map))
  ((value_type_generator)(init_type_generator))
  ((lvalue_emplacer)(norehash_lvalue_emplacer)
   (lvalue_emplace_or_cvisit)(lvalue_emplace_or_visit))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
