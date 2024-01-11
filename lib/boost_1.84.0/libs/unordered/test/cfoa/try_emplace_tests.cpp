// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

#include <boost/core/ignore_unused.hpp>

namespace {
  test::seed_t initialize_seed(511933564);

  struct lvalue_try_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto const& r : s) {
          bool b = x.try_emplace(r.first, r.second.x_);
          if (b) {
            ++num_inserts;
          }
        }
      });
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(raii::copy_constructor, x.size());
      BOOST_TEST_EQ(raii::default_constructor, x.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_try_emplacer;

  struct norehash_lvalue_try_emplacer_type : public lvalue_try_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      x.reserve(values.size());
      lvalue_try_emplacer_type::operator()(values, x);
      BOOST_TEST_EQ(raii::move_constructor, 0u);
    }
  } norehash_lvalue_try_emplacer;

  struct rvalue_try_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      BOOST_TEST_EQ(raii::copy_constructor, 0u);

      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace(std::move(r.first), r.second.x_);
          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, x.size());

      if (std::is_same<T, typename X::value_type>::value) {
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(raii::move_constructor, x.size());
      }

      BOOST_TEST_EQ(raii::default_constructor, x.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } rvalue_try_emplacer;

  struct norehash_rvalue_try_emplacer_type : public rvalue_try_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      x.reserve(values.size());

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_EQ(raii::move_constructor, 0u);

      rvalue_try_emplacer_type::operator()(values, x);

      if (std::is_same<T, typename X::value_type>::value) {
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_EQ(raii::move_constructor, x.size());
      }
    }
  } norehash_rvalue_try_emplacer;

  struct transp_try_emplace_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using is_transparent =
        typename boost::make_void<typename X::hasher::is_transparent,
          typename X::key_equal::is_transparent>::type;

      boost::ignore_unused<is_transparent>();

      BOOST_TEST_EQ(raii::default_constructor, 0u);

      std::atomic<std::uint64_t> num_inserts{0};

      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace(r.first.x_, r.second.x_);
          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(raii::default_constructor, 2 * x.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } transp_try_emplace;

  struct norehash_transp_try_emplace_type : public transp_try_emplace_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      x.reserve(values.size());
      transp_try_emplace_type::operator()(values, x);
      BOOST_TEST_EQ(raii::move_constructor, 0u);
    }
  } norehash_transp_try_emplace;

  struct lvalue_try_emplace_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace_or_cvisit(
            r.first, r.second.x_,
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

      BOOST_TEST_EQ(raii::default_constructor, x.size());
      BOOST_TEST_EQ(raii::copy_constructor, x.size());
      // don't check move construction count here because of rehashing
      BOOST_TEST_GT(raii::move_constructor, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_try_emplace_or_cvisit;

  struct lvalue_try_emplace_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace_or_visit(
            r.first, r.second.x_,
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

      BOOST_TEST_EQ(raii::default_constructor, x.size());
      BOOST_TEST_EQ(raii::copy_constructor, x.size());
      // don't check move construction count here because of rehashing
      BOOST_TEST_GT(raii::move_constructor, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_try_emplace_or_visit;

  struct rvalue_try_emplace_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace_or_cvisit(
            std::move(r.first), r.second.x_,
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

      BOOST_TEST_EQ(raii::default_constructor, x.size());

      if (std::is_same<T, typename X::value_type>::value) {
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
        BOOST_TEST_GE(raii::move_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(raii::move_constructor, x.size());
      }
    }
  } rvalue_try_emplace_or_cvisit;

  struct rvalue_try_emplace_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace_or_visit(
            std::move(r.first), r.second.x_,
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

      BOOST_TEST_EQ(raii::default_constructor, x.size());
      if (std::is_same<T, typename X::value_type>::value) {
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
        BOOST_TEST_GE(raii::move_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(raii::move_constructor, x.size());
      }
    }
  } rvalue_try_emplace_or_visit;

  struct transp_try_emplace_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace_or_cvisit(
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
      BOOST_TEST_EQ(raii::default_constructor, 2 * x.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
    }
  } transp_try_emplace_or_cvisit;

  struct transp_try_emplace_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.try_emplace_or_visit(
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

      BOOST_TEST_EQ(raii::default_constructor, 2 * x.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
    }
  } transp_try_emplace_or_visit;

  template <class X, class G, class F>
  void try_emplace(X*, G gen, F try_emplacer, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());
    raii::reset_counts();

    {
      X x;

      try_emplacer(values, x);

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
  boost::unordered::concurrent_flat_map<raii, raii, transp_hash,
    transp_key_equal>* transp_map;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

value_generator<std::pair<raii const, raii> > value_type_generator;
value_generator<std::pair<raii, raii> > init_type_generator;

// clang-format off
UNORDERED_TEST(
  try_emplace,
  ((map))
  ((value_type_generator)(init_type_generator))
  ((lvalue_try_emplacer)(norehash_lvalue_try_emplacer)
   (rvalue_try_emplacer)(norehash_rvalue_try_emplacer)
   (lvalue_try_emplace_or_cvisit)(lvalue_try_emplace_or_visit)
   (rvalue_try_emplace_or_cvisit)(rvalue_try_emplace_or_visit))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  try_emplace,
  ((transp_map))
  ((init_type_generator))
  ((transp_try_emplace)(norehash_transp_try_emplace)
   (transp_try_emplace_or_cvisit)(transp_try_emplace_or_visit))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
