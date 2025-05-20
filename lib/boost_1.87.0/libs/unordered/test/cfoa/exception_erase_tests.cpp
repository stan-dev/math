// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "exception_helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

#include <boost/core/ignore_unused.hpp>

namespace {
  test::seed_t initialize_seed(3202923);

  struct lvalue_eraser_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::atomic<std::uint64_t> num_erased{0};

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      enable_exceptions();
      thread_runner(values, [&values, &num_erased, &x](boost::span<T>) {
        for (auto const& v : values) {
          try {
            auto count = x.erase(get_key(v));
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

      BOOST_TEST_EQ(
        raii::destructor, old_d + value_type_cardinality * num_erased);
    }
  } lvalue_eraser;

  struct lvalue_eraser_if_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using value_type = typename X::value_type;
      static constexpr auto value_type_cardinality = 
        value_cardinality<value_type>::value;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      std::atomic<std::uint64_t> num_erased{0};

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      auto max = 0;
      x.visit_all([&max](value_type const& v) {
        if (get_value(v).x_ > max) {
          max = get_value(v).x_;
        }
      });

      auto threshold = max / 2;

      auto expected_erasures = 0u;
      x.visit_all([&expected_erasures, threshold](value_type const& v) {
        if (get_value(v).x_ > threshold) {
          ++expected_erasures;
        }
      });

      enable_exceptions();
      thread_runner(values, [&num_erased, &x, threshold](boost::span<T> s) {
        for (auto const& v : s) {
          try {
            auto count = x.erase_if(get_key(v),
              [threshold](arg_type& w) { return get_value(w).x_ > threshold; });
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

      BOOST_TEST_EQ(
        raii::destructor, old_d + value_type_cardinality * num_erased);
    }
  } lvalue_eraser_if;

  struct erase_if_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using value_type = typename X::value_type;
      static constexpr auto value_type_cardinality = 
        value_cardinality<value_type>::value;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      auto max = 0;
      x.visit_all([&max](value_type const& v) {
        if (get_value(v).x_ > max) {
          max = get_value(v).x_;
        }
      });

      auto threshold = max / 2;

      auto expected_erasures = 0u;
      x.visit_all([&expected_erasures, threshold](value_type const& v) {
        if (get_value(v).x_ > threshold) {
          ++expected_erasures;
        }
      });

      enable_exceptions();
      thread_runner(values, [&x, threshold](boost::span<T> /* s */) {
        for (std::size_t i = 0; i < 256; ++i) {
          try {
            x.erase_if([threshold](arg_type& v) {
              static std::atomic<std::uint32_t> c{0};
              auto t = ++c;
              if (should_throw && (t % throw_threshold == 0)) {
                throw exception_tag{};
              }

              return get_value(v).x_ > threshold;
            });
          } catch (...) {
          }
        }
      });
      disable_exceptions();

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(
        raii::destructor, 
        old_d + value_type_cardinality * (old_size - x.size()));
    }
  } erase_if;

  struct free_fn_erase_if_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using value_type = typename X::value_type;
      static constexpr auto value_type_cardinality = 
        value_cardinality<value_type>::value;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      auto const old_size = x.size();

      auto const old_dc = +raii::default_constructor;
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      auto const old_d = +raii::destructor;

      auto max = 0;
      x.visit_all([&max](value_type const& v) {
        if (get_value(v).x_ > max) {
          max = get_value(v).x_;
        }
      });

      auto threshold = max / 2;

      enable_exceptions();
      thread_runner(values, [&x, threshold](boost::span<T> /* s */) {
        for (std::size_t i = 0; i < 256; ++i) {
          try {
            boost::unordered::erase_if(x, [threshold](arg_type& v) {
              static std::atomic<std::uint32_t> c{0};
              auto t = ++c;
              if (should_throw && (t % throw_threshold == 0)) {
                throw exception_tag{};
              }

              return get_value(v).x_ > threshold;
            });

          } catch (...) {
          }
        }
      });
      disable_exceptions();

      BOOST_TEST_EQ(raii::default_constructor, old_dc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      BOOST_TEST_EQ(
        raii::destructor, 
        old_d + value_type_cardinality * (old_size - x.size()));
    }
  } free_fn_erase_if;

  template <class X, class GF, class F>
  void erase(X*, GF gen_factory, F eraser, test::random_generator rg)
  {
    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      X x(values.size());
      for (auto const& v : values) {
        x.insert(v);
      }

      BOOST_TEST_EQ(x.size(), reference_cont.size());
      BOOST_TEST_EQ(raii::destructor, 0u);

      test_fuzzy_matches_reference(x, reference_cont, rg);

      eraser(values, x);
      test_fuzzy_matches_reference(x, reference_cont, rg);
    }

    check_raii_counts();
  }

  boost::unordered::concurrent_flat_map<raii, raii, stateful_hash,
    stateful_key_equal, stateful_allocator<std::pair<raii const, raii> > >* map;
  boost::unordered::concurrent_node_map<raii, raii, stateful_hash,
    stateful_key_equal, stateful_allocator<std::pair<raii const, raii> > >* node_map;
  boost::unordered::concurrent_flat_set<raii, stateful_hash,
    stateful_key_equal, stateful_allocator<raii> >* set;
  boost::unordered::concurrent_node_set<raii, stateful_hash,
    stateful_key_equal, stateful_allocator<raii> >* node_set;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off
UNORDERED_TEST(
  erase,
  ((map)(node_map)(set)(node_set))
  ((exception_value_type_generator_factory)
   (exception_init_type_generator_factory))
  ((lvalue_eraser)(lvalue_eraser_if)(erase_if)(free_fn_erase_if))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

RUN_TESTS()
