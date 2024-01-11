// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>

#include <boost/core/ignore_unused.hpp>

namespace {
  test::seed_t initialize_seed(335740237);

  template <typename Container, typename Value>
  bool member_emplace(Container& x, Value const & v)
  {
    return x.emplace(v.x_);
  }

  template <typename Container, typename Key, typename Value>
  bool member_emplace(Container& x, std::pair<Key, Value> const & v)
  {
    return x.emplace(v.first.x_, v.second.x_);
  }

  template <typename Container, typename Value, typename F>
  bool member_emplace_or_visit(Container& x, Value& v, F f)
  {
    return x.emplace_or_visit(v.x_, f);
  }

  template <typename Container, typename Key, typename Value, typename F>
  bool member_emplace_or_visit(Container& x, std::pair<Key, Value>& v, F f)
  {
    return x.emplace_or_visit(v.first.x_, v.second.x_, f);
  }

  template <typename Container, typename Value, typename F>
  bool member_emplace_or_cvisit(Container& x, Value& v, F f)
  {
    return x.emplace_or_cvisit(v.x_, f);
  }

  template <typename Container, typename Key, typename Value, typename F>
  bool member_emplace_or_cvisit(Container& x, std::pair<Key, Value>& v, F f)
  {
    return x.emplace_or_cvisit(v.first.x_, v.second.x_, f);
  }

  struct lvalue_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto const& r : s) {
          bool b = member_emplace(x, r);
          if (b) {
            ++num_inserts;
          }
        }
      });
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values.size());

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, value_type_cardinality * x.size());

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_emplacer;

  struct norehash_lvalue_emplacer_type : public lvalue_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      x.reserve(values.size());
      lvalue_emplacer_type::operator()(values, x);
      BOOST_TEST_EQ(raii::move_constructor, value_type_cardinality * x.size());
    }
  } norehash_lvalue_emplacer;

  struct lvalue_emplace_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = member_emplace_or_cvisit(
            x, r,
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

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_emplace_or_cvisit;

  struct lvalue_emplace_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = member_emplace_or_visit(
            x, r,
            [&num_invokes](arg_type& v) {
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

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_emplace_or_visit;

  template <class X, class GF, class F>
  void emplace(X*, GF gen_factory, F emplacer, test::random_generator rg)
  {
    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());
    raii::reset_counts();

    {
      X x;

      emplacer(values, x);

      BOOST_TEST_EQ(x.size(), reference_cont.size());

      using value_type = typename X::value_type;
      BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& v) {
        BOOST_TEST(reference_cont.contains(get_key(v)));
        if (rg == test::sequential) {
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
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
  boost::unordered::concurrent_flat_set<raii>* set;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off

UNORDERED_TEST(
  emplace,
  ((map)(set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((lvalue_emplacer)(norehash_lvalue_emplacer)
   (lvalue_emplace_or_cvisit)(lvalue_emplace_or_visit))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
