// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/config.hpp>
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

#include <boost/core/ignore_unused.hpp>

#if defined(BOOST_MSVC)
#pragma warning(disable : 4127) // conditional expression is constant
#endif

struct raii_convertible
{
  int x = 0, y = 0 ;

  template <typename T>
  raii_convertible(T const & t) : x{t.x_} {}

  template <typename T, typename Q>
  raii_convertible(std::pair<T, Q> const & p) : x{p.first.x_}, y{p.second.x_} 
  {}

  operator raii() { return {x}; }
  operator std::pair<raii const, raii>() { return {x, y}; }
};

namespace {
  test::seed_t initialize_seed(78937);

  struct lvalue_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto const& r : s) {
          bool b = x.insert(r);
          if (b) {
            ++num_inserts;
          }
        }
      });
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_inserter;

  struct norehash_lvalue_inserter_type : public lvalue_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      x.reserve(values.size());
      lvalue_inserter_type::operator()(values, x);
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::move_constructor, 0u);
    }
  } norehash_lvalue_inserter;

  struct rvalue_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      BOOST_TEST_EQ(raii::copy_constructor, 0u);

      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.insert(std::move(r));
          if (b) {
            ++num_inserts;
          }
        }
      });
      BOOST_TEST_EQ(num_inserts, x.size());

      if (std::is_same<T, typename X::value_type>::value &&
          !std::is_same<typename X::key_type, typename X::value_type>::value) {
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
      }

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } rvalue_inserter;

  struct norehash_rvalue_inserter_type : public rvalue_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      x.reserve(values.size());

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_EQ(raii::move_constructor, 0u);

      rvalue_inserter_type::operator()(values, x);

      if (std::is_same<T, typename X::value_type>::value) {
        if (std::is_same<typename X::key_type, 
                         typename X::value_type>::value) {
          BOOST_TEST_EQ(raii::copy_constructor, 0u);
          BOOST_TEST_EQ(raii::move_constructor, x.size());
        }
        else {
          BOOST_TEST_EQ(raii::copy_constructor, x.size());
          BOOST_TEST_EQ(raii::move_constructor, x.size());
        }
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_EQ(
          raii::move_constructor, value_type_cardinality * x.size());
      }
    }
  } norehash_rvalue_inserter;

  struct iterator_range_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& v : values) { 
        values2.push_back(raii_convertible(v));
      }

      thread_runner(values2, [&x](boost::span<raii_convertible> s) {
        BOOST_TEST_EQ(x.insert(s.begin(), s.end()), s.size());
      });

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values2.size());
#if BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
    BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)
      // some versions of old gcc have trouble eliding copies here
      // https://godbolt.org/z/Ebo6TbvaG
#elif BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 40900) && \
      BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50000)
      // seemingly same problem, though the snippet above does not reveal it
#else
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
#endif
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } iterator_range_inserter;

  struct lvalue_insert_or_assign_copy_assign_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      thread_runner(values, [&x](boost::span<T> s) {
        for (auto& r : s) {
          x.insert_or_assign(r.first, r.second);
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 2 * x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      }
      else{
        // don't check move construction count here because of rehashing
        BOOST_TEST_GT(raii::move_constructor, 0u);
      }

      BOOST_TEST_EQ(raii::copy_assignment, values.size() - x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_insert_or_assign_copy_assign;

  struct lvalue_insert_or_assign_move_assign_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      thread_runner(values, [&x](boost::span<T> s) {
        for (auto& r : s) {
          x.insert_or_assign(r.first, std::move(r.second));
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, x.size());
      }
      else{
        BOOST_TEST_GT(raii::move_constructor, x.size()); // rehashing
      }

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, values.size() - x.size());
    }
  } lvalue_insert_or_assign_move_assign;

  struct rvalue_insert_or_assign_copy_assign_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      thread_runner(values, [&x](boost::span<T> s) {
        for (auto& r : s) {
          x.insert_or_assign(std::move(r.first), r.second);
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, x.size());
      }
      else{
        BOOST_TEST_GT(raii::move_constructor, x.size()); // rehashing
      }

      BOOST_TEST_EQ(raii::copy_assignment, values.size() - x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } rvalue_insert_or_assign_copy_assign;

  struct rvalue_insert_or_assign_move_assign_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      thread_runner(values, [&x](boost::span<T> s) {
        for (auto& r : s) {
          x.insert_or_assign(std::move(r.first), std::move(r.second));
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 0u);

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 2 * x.size());
      }
      else{
        BOOST_TEST_GE(raii::move_constructor, 2 * x.size()); // rehashing
      }

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, values.size() - x.size());
    }
  } rvalue_insert_or_assign_move_assign;

  struct trans_insert_or_assign_copy_assign_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using is_transparent =
        typename boost::make_void<typename X::hasher::is_transparent,
          typename X::key_equal::is_transparent>::type;

      boost::ignore_unused<is_transparent>();

      BOOST_TEST_EQ(raii::default_constructor, 0u);

      thread_runner(values, [&x](boost::span<T> s) {
        for (auto& r : s) {
          x.insert_or_assign(r.first.x_, r.second);
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, x.size());
      BOOST_TEST_EQ(raii::copy_constructor, x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      }
      else{
        BOOST_TEST_GT(raii::move_constructor, 0u); // rehashing
      }

      BOOST_TEST_EQ(raii::copy_assignment, values.size() - x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } trans_insert_or_assign_copy_assign;

  struct trans_insert_or_assign_move_assign_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      using is_transparent =
        typename boost::make_void<typename X::hasher::is_transparent,
          typename X::key_equal::is_transparent>::type;

      boost::ignore_unused<is_transparent>();

      thread_runner(values, [&x](boost::span<T> s) {
        for (auto& r : s) {
          x.insert_or_assign(r.first.x_, std::move(r.second));
        }
      });

      BOOST_TEST_EQ(raii::default_constructor, x.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, x.size());
      }
      else{
        BOOST_TEST_GT(raii::move_constructor, x.size()); // rehashing
      }

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, values.size() - x.size());
    }
  } trans_insert_or_assign_move_assign;

  struct lvalue_insert_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.insert_or_cvisit(
            r, [&num_invokes](typename X::value_type const& v) {
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

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      }
      else{
        // don't check move construction count here because of rehashing
        BOOST_TEST_GT(raii::move_constructor, 0u);
      }

      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_insert_or_cvisit;

  struct lvalue_insert_and_cvisit_type
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

      std::atomic<std::uint64_t> num_inserts{0}, num_inserts_internal{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values,
        [&x, &num_inserts, &num_inserts_internal, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.insert_and_cvisit(
            r,
            [&num_inserts_internal](arg_type& v) {
              (void)v;
              ++num_inserts_internal;
            },
            [&num_invokes](typename X::value_type const& v) {
              (void)v;
              ++num_invokes;
            });

          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, num_inserts_internal);
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      }
      else{
        // don't check move construction count here because of rehashing
        BOOST_TEST_GT(raii::move_constructor, 0u);
      }

      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_insert_and_cvisit;

  struct lvalue_insert_or_visit_type
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
          bool b =
            x.insert_or_visit(r, [&num_invokes](arg_type& v) {
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

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, value_type_cardinality * x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      }
      else{
        // don't check move construction count here because of rehashing
        BOOST_TEST_GT(raii::move_constructor, 0u);
      }

      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_insert_or_visit;

  struct lvalue_insert_and_visit_type
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

      std::atomic<std::uint64_t> num_inserts{0}, num_inserts_internal{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values,
        [&x, &num_inserts, &num_inserts_internal, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b =
            x.insert_and_visit(r, 
            [&num_inserts_internal](arg_type& v) {
              (void)v;
              ++num_inserts_internal;
            },
            [&num_invokes](arg_type& v) {
              (void)v;
              ++num_invokes;
            });

          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, num_inserts_internal);
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, value_type_cardinality * x.size());

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      }
      else{
        // don't check move construction count here because of rehashing
        BOOST_TEST_GT(raii::move_constructor, 0u);
      }

      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_insert_and_visit;

  struct rvalue_insert_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.insert_or_cvisit(
            std::move(r), [&num_invokes](typename X::value_type const& v) {
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

      BOOST_TEST_EQ(raii::default_constructor, 0u);

      if (std::is_same<T, typename X::value_type>::value) {
        if (std::is_same<typename X::key_type, 
                         typename X::value_type>::value) {
          BOOST_TEST_EQ(raii::copy_constructor, 0u);
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
        else {
          BOOST_TEST_EQ(raii::copy_constructor, x.size());
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(
          raii::move_constructor, value_type_cardinality * x.size());
      }
    }
  } rvalue_insert_or_cvisit;

  struct rvalue_insert_and_cvisit_type
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

      std::atomic<std::uint64_t> num_inserts{0}, num_inserts_internal{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, 
        [&x, &num_inserts, &num_inserts_internal, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.insert_and_cvisit(
            std::move(r), 
            [&num_inserts_internal](arg_type& v) {
              (void)v;
              ++num_inserts_internal;
            },            
            [&num_invokes](typename X::value_type const& v) {
              (void)v;
              ++num_invokes;
            });

          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, num_inserts_internal);
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 0u);

      if (std::is_same<T, typename X::value_type>::value) {
        if (std::is_same<typename X::key_type, 
                         typename X::value_type>::value) {
          BOOST_TEST_EQ(raii::copy_constructor, 0u);
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
        else {
          BOOST_TEST_EQ(raii::copy_constructor, x.size());
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(
          raii::move_constructor, value_type_cardinality * x.size());
      }
    }
  } rvalue_insert_and_cvisit;

  struct rvalue_insert_or_visit_type
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
          bool b = x.insert_or_visit(
            std::move(r), [&num_invokes](arg_type& v) {
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

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      if (std::is_same<T, typename X::value_type>::value) {
        if (std::is_same<typename X::key_type, 
                         typename X::value_type>::value) {
          BOOST_TEST_EQ(raii::copy_constructor, 0u);
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
        else {
          BOOST_TEST_EQ(raii::copy_constructor, x.size());
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(
          raii::move_constructor, value_type_cardinality * x.size());
      }
    }
  } rvalue_insert_or_visit;

  struct rvalue_insert_and_visit_type
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

      std::atomic<std::uint64_t> num_inserts{0}, num_inserts_internal{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values,
        [&x, &num_inserts, &num_inserts_internal, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.insert_and_visit(
            std::move(r), 
            [&num_inserts_internal](arg_type& v) {
              (void)v;
              ++num_inserts_internal;
            },            
            [&num_invokes](arg_type& v) {
              (void)v;
              ++num_invokes;
            });

          if (b) {
            ++num_inserts;
          }
        }
      });

      BOOST_TEST_EQ(num_inserts, num_inserts_internal);
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      if (std::is_same<T, typename X::value_type>::value) {
        if (std::is_same<typename X::key_type, 
                         typename X::value_type>::value) {
          BOOST_TEST_EQ(raii::copy_constructor, 0u);
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
        else {
          BOOST_TEST_EQ(raii::copy_constructor, x.size());
          BOOST_TEST_GE(raii::move_constructor, x.size());
        }
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(
          raii::move_constructor, value_type_cardinality * x.size());
      }
    }
  } rvalue_insert_and_visit;

  struct iterator_range_insert_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality = 
        value_cardinality<typename X::value_type>::value;

      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& v : values) {
        values2.push_back(raii_convertible(v));
      }

      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(
        values2, [&x, &num_invokes](boost::span<raii_convertible> s) {
          BOOST_TEST_EQ(x.insert_or_cvisit(s.begin(), s.end(),
                          [&num_invokes](typename X::value_type const& v) {
                            (void)v;
                            ++num_invokes;
                          }),
            s.size());
        });

      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values2.size());
#if (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)) || \
    (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 40900) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50000))
      // skip test
#else
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
#endif
      BOOST_TEST_GT(raii::move_constructor, 0u);
    }
  } iterator_range_insert_or_cvisit;

  struct iterator_range_insert_and_cvisit_type
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

      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& v : values) {
        values2.push_back(raii_convertible(v));
      }

      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values2,
        [&x, &num_inserts, &num_invokes](boost::span<raii_convertible> s) {
          BOOST_TEST_EQ(x.insert_and_cvisit(
                          s.begin(), s.end(),
                          [&num_inserts](arg_type& v) {
                            (void)v;
                            ++num_inserts;
                          },
                          [&num_invokes](typename X::value_type const& v) {
                            (void)v;
                            ++num_invokes;
                          }),
            s.size());
        });

      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values2.size());
#if (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)) || \
    (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 40900) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50000))
      // skip test
#else
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
#endif
      BOOST_TEST_GT(raii::move_constructor, 0u);
    }
  } iterator_range_insert_and_cvisit;

  struct iterator_range_insert_or_visit_type
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

      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& v : values) {
        values2.push_back(raii_convertible(v));
      }

      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(
        values2, [&x, &num_invokes](boost::span<raii_convertible> s) {
          BOOST_TEST_EQ(x.insert_or_visit(s.begin(), s.end(),
                          [&num_invokes](arg_type& v) {
                            (void)v;
                            ++num_invokes;
                          }),
            s.size());
        });

      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values2.size());
#if (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)) || \
    (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 40900) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50000))
      // skip test
#else
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
#endif
      BOOST_TEST_GT(raii::move_constructor, 0u);
    }
  } iterator_range_insert_or_visit;

  struct iterator_range_insert_and_visit_type
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

      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& v : values) {
        values2.push_back(raii_convertible(v));
      }

      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values2,
        [&x, &num_inserts, &num_invokes](boost::span<raii_convertible> s) {
          BOOST_TEST_EQ(x.insert_and_visit(
                          s.begin(), s.end(),
                          [&num_inserts](arg_type& v) {
                            (void)v;
                            ++num_inserts;
                          },
                          [&num_invokes](typename X::value_type const& v) {
                            (void)v;
                            ++num_invokes;
                          }),
            s.size());
        });

      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values2.size());
#if (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)) || \
    (BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 40900) && \
     BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50000))
      // skip test
#else
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
#endif
      BOOST_TEST_GT(raii::move_constructor, 0u);
    }
  } iterator_range_insert_and_visit;

  struct non_copyable_function
  {
    non_copyable_function() = default;
    non_copyable_function(const non_copyable_function&) = delete;
    non_copyable_function(non_copyable_function&&) = default;
    template <class... Args> void operator()(Args&&...) const {}
  };

  template <class X, class GF, class F>
  void insert(X*, GF gen_factory, F inserter, test::random_generator rg)
  {
    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());
    raii::reset_counts();

    {
      X x;

      inserter(values, x);

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

  template <class X, class IL>
  void insert_initializer_list(std::pair<X*, IL> p)
  {
    using value_type = typename X::value_type;
    // concurrent_flat_set visit is always const access
    using arg_type = typename std::conditional<
      std::is_same<typename X::key_type, typename X::value_type>::value,
      typename X::value_type const,
      typename X::value_type
    >::type;

    auto init_list = p.second;
    std::vector<raii> dummy;
    auto reference_cont = reference_container<X>(
      init_list.begin(), init_list.end());
    raii::reset_counts();

    {
      {
        X x;

        thread_runner(dummy, [&x, &init_list](boost::span<raii>) {
          BOOST_TEST_EQ(x.insert(init_list), init_list.size());
        });

        BOOST_TEST_EQ(x.size(), reference_cont.size());

        BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
        }));
      }

      BOOST_TEST_GE(raii::default_constructor, 0u);
      BOOST_TEST_GE(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, 0u);
      BOOST_TEST_GT(raii::destructor, 0u);

      BOOST_TEST_EQ(raii::default_constructor + raii::copy_constructor +
                      raii::move_constructor,
        raii::destructor);

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }

    {
      {
        std::atomic<std::uint64_t> num_invokes{0};

        X x;

        thread_runner(dummy, [&x, &init_list, &num_invokes](boost::span<raii>) {
          BOOST_TEST_EQ(x.insert_or_visit(init_list,
                          [&num_invokes](arg_type& v) {
                            (void)v;
                            ++num_invokes;
                          }),
            init_list.size());

          BOOST_TEST_EQ(x.insert_or_cvisit(init_list,
                          [&num_invokes](typename X::value_type const& v) {
                            (void)v;
                            ++num_invokes;
                          }),
            init_list.size());

          x.insert_or_visit(init_list, non_copyable_function{});
          x.insert_or_cvisit(init_list, non_copyable_function{});
        });

        BOOST_TEST_EQ(num_invokes, (init_list.size() - x.size()) +
                                     (num_threads - 1) * init_list.size() +
                                     num_threads * init_list.size());
        BOOST_TEST_EQ(x.size(), reference_cont.size());

        BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
        }));
      }

      BOOST_TEST_GE(raii::default_constructor, 0u);
      BOOST_TEST_GE(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, 0u);
      BOOST_TEST_GT(raii::destructor, 0u);

      BOOST_TEST_EQ(raii::default_constructor + raii::copy_constructor +
                      raii::move_constructor,
        raii::destructor);

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }

    {
      {
        std::atomic<std::uint64_t> num_inserts{0};
        std::atomic<std::uint64_t> num_invokes{0};

        X x;

        thread_runner(dummy,
          [&x, &init_list, &num_inserts, &num_invokes](boost::span<raii>) {
            BOOST_TEST_EQ(x.insert_and_visit(
                            init_list,
                            [&num_inserts](arg_type& v) {
                              (void)v;
                              ++num_inserts;
                            },
                            [&num_invokes](arg_type& v) {
                              (void)v;
                              ++num_invokes;
                            }),
              init_list.size());

            BOOST_TEST_EQ(x.insert_and_cvisit(
                            init_list,
                            [&num_inserts](arg_type& v) {
                              (void)v;
                              ++num_inserts;
                            },
                            [&num_invokes](typename X::value_type const& v) {
                              (void)v;
                              ++num_invokes;
                            }),
              init_list.size());

            x.insert_and_visit(
              init_list, non_copyable_function{}, non_copyable_function{});
            x.insert_and_cvisit(
              init_list, non_copyable_function{}, non_copyable_function{});
          });

        BOOST_TEST_EQ(num_inserts, x.size());
        BOOST_TEST_EQ(num_invokes, (init_list.size() - x.size()) +
                                     (num_threads - 1) * init_list.size() +
                                     num_threads * init_list.size());
        BOOST_TEST_EQ(x.size(), reference_cont.size());

        BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
        }));
      }

      BOOST_TEST_GE(raii::default_constructor, 0u);
      BOOST_TEST_GE(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, 0u);
      BOOST_TEST_GT(raii::destructor, 0u);

      BOOST_TEST_EQ(raii::default_constructor + raii::copy_constructor +
                      raii::move_constructor,
        raii::destructor);

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  }


  template <class X>
  void insert_map_sfinae_test(X*)
  {
    // mostly a compile-time tests to ensure that there's no ambiguity when a
    // user does this
    using value_type = typename X::value_type;
    X x;
    x.insert({1, 2});

    x.insert_or_visit({2, 3}, [](value_type&) {});
    x.insert_or_cvisit({3, 4}, [](value_type const&) {});

    x.insert_and_visit({4, 5}, [](value_type&) {}, [](value_type&) {});
    x.insert_and_cvisit({5, 6}, [](value_type&) {}, [](value_type const&) {});
  }

  boost::unordered::concurrent_flat_map<raii, raii>* map;
  boost::unordered::concurrent_flat_map<raii, raii, transp_hash,
    transp_key_equal>* trans_map;
  boost::unordered::concurrent_flat_map<raii, raii, boost::hash<raii>,
    std::equal_to<raii>, fancy_allocator<std::pair<raii const, raii> > >*
    fancy_map;

  boost::unordered::concurrent_node_map<raii, raii>* node_map;
  boost::unordered::concurrent_node_map<raii, raii, transp_hash,
    transp_key_equal>* trans_node_map;
  boost::unordered::concurrent_node_map<raii, raii, boost::hash<raii>,
    std::equal_to<raii>, fancy_allocator<std::pair<raii const, raii> > >*
    fancy_node_map;

  boost::unordered::concurrent_flat_set<raii>* set;
  boost::unordered::concurrent_flat_set<raii, boost::hash<raii>,
    std::equal_to<raii>, fancy_allocator<raii> >*
    fancy_set;

  boost::unordered::concurrent_node_set<raii>* node_set;
  boost::unordered::concurrent_node_set<raii, boost::hash<raii>,
    std::equal_to<raii>, fancy_allocator<raii> >*
    fancy_node_set;

  std::initializer_list<std::pair<raii const, raii> > map_init_list{
    {raii{0}, raii{0}},
    {raii{1}, raii{1}},
    {raii{2}, raii{2}},
    {raii{3}, raii{3}},
    {raii{4}, raii{4}},
    {raii{5}, raii{5}},
    {raii{6}, raii{6}},
    {raii{6}, raii{6}},
    {raii{7}, raii{7}},
    {raii{8}, raii{8}},
    {raii{9}, raii{9}},
    {raii{10}, raii{10}},
    {raii{9}, raii{9}},
    {raii{8}, raii{8}},
    {raii{7}, raii{7}},
    {raii{6}, raii{6}},
    {raii{5}, raii{5}},
    {raii{4}, raii{4}},
    {raii{3}, raii{3}},
    {raii{2}, raii{2}},
    {raii{1}, raii{1}},
    {raii{0}, raii{0}},
  };

  std::initializer_list<raii> set_init_list{
    raii{0},
    raii{1},
    raii{2},
    raii{3},
    raii{4},
    raii{5},
    raii{6},
    raii{6},
    raii{7},
    raii{8},
    raii{9},
    raii{10},
    raii{9},
    raii{8},
    raii{7},
    raii{6},
    raii{5},
    raii{4},
    raii{3},
    raii{2},
    raii{1},
    raii{0},
  };

  auto map_and_init_list=std::make_pair(map,map_init_list);  
  auto node_map_and_init_list=std::make_pair(node_map,map_init_list);  
  auto set_and_init_list=std::make_pair(set,set_init_list);
  auto node_set_and_init_list=std::make_pair(node_set,set_init_list);

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off
UNORDERED_TEST(
  insert_initializer_list,
  ((map_and_init_list)(node_map_and_init_list)
   (set_and_init_list)(node_set_and_init_list)))

UNORDERED_TEST(
  insert,
  ((map)(fancy_map)(node_map)(fancy_node_map)
   (set)(fancy_set)(node_set)(fancy_node_set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((lvalue_inserter)(rvalue_inserter)(iterator_range_inserter)
   (norehash_lvalue_inserter)(norehash_rvalue_inserter)
   (lvalue_insert_or_cvisit)(lvalue_insert_or_visit)
   (rvalue_insert_or_cvisit)(rvalue_insert_or_visit)
   (iterator_range_insert_or_cvisit)(iterator_range_insert_or_visit))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  insert,
  ((map)(fancy_map)(node_map)(fancy_node_map)
   (set)(fancy_set)(node_set)(fancy_node_set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((lvalue_insert_and_cvisit)(lvalue_insert_and_visit)
   (rvalue_insert_and_cvisit)(rvalue_insert_and_visit)
   (iterator_range_insert_and_cvisit)(iterator_range_insert_and_visit))
  ((default_generator)(sequential)(limited_range)))

 UNORDERED_TEST(
  insert,
  ((map)(node_map))
  ((init_type_generator_factory))
  ((lvalue_insert_or_assign_copy_assign)(lvalue_insert_or_assign_move_assign)
   (rvalue_insert_or_assign_copy_assign)(rvalue_insert_or_assign_move_assign))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  insert,
  ((trans_map)(trans_node_map))
  ((init_type_generator_factory))
  ((trans_insert_or_assign_copy_assign)(trans_insert_or_assign_move_assign))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  insert_map_sfinae_test,
  ((map)(node_map)))
// clang-format on

RUN_TESTS()
