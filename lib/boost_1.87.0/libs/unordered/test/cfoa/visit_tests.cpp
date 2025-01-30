// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>
#include <boost/config/workaround.hpp>

#if BOOST_WORKAROUND(BOOST_GCC_VERSION, < 40900)
// warning triggered in transform_iterator.hpp transitive includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

#include <boost/compat/latch.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/iterator/transform_iterator.hpp>

#if BOOST_WORKAROUND(BOOST_GCC_VERSION, < 40900)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <iterator>
#include <thread>
#include <vector>

namespace {
  test::seed_t initialize_seed(335740237);

  auto non_present_keys = []
  {
    std::array<raii,128> a;
    for(std::size_t i = 0; i < a.size(); ++i) {
      a[i].x_ = -((int)i + 1);
    }
    return a;
  }();

  template<typename T>
  raii const & get_non_present_key(T const & x)
  {
    return non_present_keys[
      (std::size_t)get_key(x).x_ % non_present_keys.size()];
  }

  struct lvalue_visitor_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_cont)
    {
      using value_type = typename X::value_type;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      std::atomic<std::uint64_t> num_visits{0};
      std::atomic<std::uint64_t> total_count{0};

      auto mut_visitor = [&num_visits, &reference_cont](arg_type& v) {
        BOOST_TEST(reference_cont.contains(get_key(v)));
        BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
        ++num_visits;
      };

      auto const_visitor =
        [&num_visits, &reference_cont](value_type const& v) {
        BOOST_TEST(reference_cont.contains(get_key(v)));
        BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
        ++num_visits;
      };

      {
        thread_runner(
          values, [&x, &mut_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = get_key(val).x_;
              BOOST_TEST(r >= 0);

              auto count = x.visit(get_key(val), mut_visitor);
              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.visit(get_non_present_key(val), mut_visitor);
              BOOST_TEST_EQ(count, 0u);
            }
          });

        BOOST_TEST_EQ(num_visits, values.size());
        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(
          values, [&x, &const_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = get_key(val).x_;
              BOOST_TEST(r >= 0);

              auto const& y = x;
              auto count = y.visit(get_key(val), const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = y.visit(get_non_present_key(val), const_visitor);
              BOOST_TEST_EQ(count, 0u);
            }
          });

        BOOST_TEST_EQ(num_visits, values.size());
        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(
          values, [&x, &const_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = get_key(val).x_;
              BOOST_TEST(r >= 0);

              auto count = x.cvisit(get_key(val), const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.cvisit(get_non_present_key(val), const_visitor);
              BOOST_TEST_EQ(count, 0u);
            }
          });

        BOOST_TEST_EQ(num_visits, values.size());
        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(values, [&x, &total_count](boost::span<T> s) {
          for (auto const& val : s) {
            auto r = get_key(val).x_;
            BOOST_TEST(r >= 0);

            auto count = x.count(get_key(val));
            BOOST_TEST_EQ(count, 1u);
            total_count += count;

            count = x.count(get_non_present_key(val));
            BOOST_TEST_EQ(count, 0u);
          }
        });

        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(values, [&x](boost::span<T> s) {
          for (auto const& val : s) {
            auto r = get_key(val).x_;
            BOOST_TEST(r >= 0);

            auto contains = x.contains(get_key(val));
            BOOST_TEST(contains);

            contains = x.contains(get_non_present_key(val));
            BOOST_TEST(!contains);
          }
        });

        num_visits = 0;
        total_count = 0;
      }
    }
  } lvalue_visitor;

  struct transp_visitor_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_cont)
    {
      using value_type = typename X::value_type;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      std::atomic<std::uint64_t> num_visits{0};
      std::atomic<std::uint64_t> total_count{0};

      auto mut_visitor = [&num_visits, &reference_cont](arg_type& v) {
        BOOST_TEST(reference_cont.contains(get_key(v)));
        BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
        ++num_visits;
      };

      auto const_visitor = [&num_visits, &reference_cont](value_type const& v) {
        BOOST_TEST(reference_cont.contains(get_key(v)));
        BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
        ++num_visits;
      };

      {
        thread_runner(
          values, [&x, &mut_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = get_key(val).x_;
              BOOST_TEST(r >= 0);

              auto count = x.visit(get_key(val).x_, mut_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.visit(get_non_present_key(val).x_, mut_visitor);
              BOOST_TEST_EQ(count, 0u);
            }
          });

        BOOST_TEST_EQ(num_visits, values.size());
        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(
          values, [&x, &const_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = get_key(val).x_;
              BOOST_TEST(r >= 0);

              auto const& y = x;
              auto count = y.visit(get_key(val).x_, const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = y.visit(get_non_present_key(val).x_, const_visitor);
              BOOST_TEST_EQ(count, 0u);
            }
          });

        BOOST_TEST_EQ(num_visits, values.size());
        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(
          values, [&x, &const_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = get_key(val).x_;
              BOOST_TEST(r >= 0);

              auto count = x.cvisit(get_key(val).x_, const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.cvisit(get_non_present_key(val).x_, const_visitor);
              BOOST_TEST_EQ(count, 0u);
            }
          });

        BOOST_TEST_EQ(num_visits, values.size());
        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(values, [&x, &total_count](boost::span<T> s) {
          for (auto const& val : s) {
            auto r = get_key(val).x_;
            BOOST_TEST(r >= 0);

            auto count = x.count(get_key(val).x_);
            BOOST_TEST_EQ(count, 1u);
            total_count += count;

            count = x.count(get_non_present_key(val).x_);
            BOOST_TEST_EQ(count, 0u);
          }
        });

        BOOST_TEST_EQ(total_count, values.size());

        num_visits = 0;
        total_count = 0;
      }

      {
        thread_runner(values, [&x](boost::span<T> s) {
          for (auto const& val : s) {
            auto r = get_key(val).x_;
            BOOST_TEST(r >= 0);

            auto contains = x.contains(get_key(val).x_);
            BOOST_TEST(contains);

            contains = x.contains(get_non_present_key(val).x_);
            BOOST_TEST(!contains);
          }
        });

        num_visits = 0;
        total_count = 0;
      }
    }
  } transp_visitor;

  struct visit_all_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_cont)
    {
      using value_type = typename X::value_type;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      std::atomic<std::uint64_t> total_count{0};

      auto mut_visitor = [&reference_cont](std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](arg_type& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
        };
      };

      auto const_visitor = [&reference_cont](std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
        };
      };

      {
        thread_runner(values, [&x, &total_count, &mut_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          total_count += x.visit_all(mut_visitor(num_visits));
          BOOST_TEST_EQ(x.size(), num_visits);
        });

        BOOST_TEST_EQ(total_count, num_threads * x.size());
        total_count = 0;
      }

      {
        thread_runner(
          values, [&x, &total_count, &const_visitor](boost::span<T>) {
            std::atomic<std::uint64_t> num_visits{0};
            auto const& y = x;
            total_count += y.visit_all(const_visitor(num_visits));
            BOOST_TEST_EQ(x.size(), num_visits);
          });

        BOOST_TEST_EQ(total_count, num_threads * x.size());
        total_count = 0;
      }

      {
        thread_runner(
          values, [&x, &total_count, &const_visitor](boost::span<T>) {
            std::atomic<std::uint64_t> num_visits{0};
            total_count += x.cvisit_all(const_visitor(num_visits));
            BOOST_TEST_EQ(x.size(), num_visits);
          });

        BOOST_TEST_EQ(total_count, num_threads * x.size());
        total_count = 0;
      }
    }

  } visit_all;

  struct visit_while_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_cont)
    {
      using value_type = typename X::value_type;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      auto mut_truthy_visitor = [&reference_cont](
                                  std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](arg_type& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
          return true;
        };
      };

      auto const_truthy_visitor = [&reference_cont](
                                    std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
          return true;
        };
      };

      auto mut_falsey_visitor = [&reference_cont](
                                  std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](arg_type& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          ++num_visits;
          return (get_value(v).x_ % 100) == 0;
        };
      };

      auto const_falsey_visitor = [&reference_cont](
                                    std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          ++num_visits;
          return (get_value(v).x_ % 100) == 0;
        };
      };

      {
        thread_runner(values, [&x, &mut_truthy_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST(x.visit_while(mut_truthy_visitor(num_visits)));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &const_truthy_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          auto const& y = x;
          BOOST_TEST(y.visit_while(const_truthy_visitor(num_visits)));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &const_truthy_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST(x.cvisit_while(const_truthy_visitor(num_visits)));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &mut_falsey_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST_NOT(x.visit_while(mut_falsey_visitor(num_visits)));
          BOOST_TEST_LT(num_visits, x.size());
          BOOST_TEST_GT(num_visits, 0u);
        });
      }

      {
        thread_runner(values, [&x, &const_falsey_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          auto const& y = x;
          BOOST_TEST_NOT(y.visit_while(const_falsey_visitor(num_visits)));
          BOOST_TEST_LT(num_visits, x.size());
          BOOST_TEST_GT(num_visits, 0u);
        });
      }

      {
        thread_runner(values, [&x, &const_falsey_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST_NOT(x.cvisit_while(const_falsey_visitor(num_visits)));
          BOOST_TEST_LT(num_visits, x.size());
          BOOST_TEST_GT(num_visits, 0u);
        });
      }
    }
  } visit_while;

  struct exec_policy_visit_all_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_cont)
    {
#if defined(BOOST_UNORDERED_PARALLEL_ALGORITHMS)
      using value_type = typename X::value_type;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      auto mut_visitor = [&reference_cont](std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](arg_type& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
        };
      };

      auto const_visitor = [&reference_cont](std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
        };
      };

      {
        thread_runner(values, [&x, &mut_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};

          x.visit_all(std::execution::par, mut_visitor(num_visits));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &const_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          auto const& y = x;

          y.visit_all(std::execution::par, const_visitor(num_visits));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &const_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          x.cvisit_all(std::execution::par, const_visitor(num_visits));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }
#else
      (void)values;
      (void)x;
      (void)reference_cont;
#endif
    }
  } exec_policy_visit_all;

  struct exec_policy_visit_while_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_cont)
    {
#if defined(BOOST_UNORDERED_PARALLEL_ALGORITHMS)
      using value_type = typename X::value_type;

      // concurrent_flat_set visit is always const access
      using arg_type = typename std::conditional<
        std::is_same<typename X::key_type, typename X::value_type>::value,
        typename X::value_type const,
        typename X::value_type
      >::type;

      auto mut_truthy_visitor = [&reference_cont](
                                  std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](arg_type& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
          return true;
        };
      };

      auto const_truthy_visitor = [&reference_cont](
                                    std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
          return true;
        };
      };

      auto mut_falsey_visitor = [&reference_cont](
                                  std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](arg_type& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
          return (get_value(v).x_ % 100) == 0;
        };
      };

      auto const_falsey_visitor = [&reference_cont](
                                    std::atomic<uint64_t>& num_visits) {
        return [&reference_cont, &num_visits](value_type const& v) {
          BOOST_TEST(reference_cont.contains(get_key(v)));
          BOOST_TEST_EQ(v, *reference_cont.find(get_key(v)));
          ++num_visits;
          return (get_value(v).x_ % 100) == 0;
        };
      };

      {
        thread_runner(values, [&x, &mut_truthy_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST(
            x.visit_while(std::execution::par, mut_truthy_visitor(num_visits)));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &const_truthy_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          auto const& y = x;
          BOOST_TEST(y.visit_while(
            std::execution::par, const_truthy_visitor(num_visits)));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &const_truthy_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST(x.cvisit_while(
            std::execution::par, const_truthy_visitor(num_visits)));
          BOOST_TEST_EQ(x.size(), num_visits);
        });
      }

      {
        thread_runner(values, [&x, &mut_falsey_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST_NOT(
            x.visit_while(std::execution::par, mut_falsey_visitor(num_visits)));
          BOOST_TEST_LT(num_visits, x.size());
          BOOST_TEST_GT(num_visits, 0u);
        });
      }

      {
        thread_runner(values, [&x, &const_falsey_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          auto const& y = x;
          BOOST_TEST_NOT(y.visit_while(
            std::execution::par, const_falsey_visitor(num_visits)));
          BOOST_TEST_LT(num_visits, x.size());
          BOOST_TEST_GT(num_visits, 0u);
        });
      }

      {
        thread_runner(values, [&x, &const_falsey_visitor](boost::span<T>) {
          std::atomic<std::uint64_t> num_visits{0};
          BOOST_TEST_NOT(x.cvisit_while(
            std::execution::par, const_falsey_visitor(num_visits)));
          BOOST_TEST_LT(num_visits, x.size());
          BOOST_TEST_GT(num_visits, 0u);
        });
      }
#else
      (void)values;
      (void)x;
      (void)reference_cont;
#endif
    }
  } exec_policy_visit_while;

  template <class X, class GF, class F>
  void visit(X*, GF gen_factory, F visitor, test::random_generator rg)
  {
    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      X x;
      for (auto const& v : values) {
        x.insert(v);
      }
      BOOST_TEST_EQ(x.size(), reference_cont.size());

      std::uint64_t old_default_constructor = raii::default_constructor;
      std::uint64_t old_copy_constructor = raii::copy_constructor;
      std::uint64_t old_move_constructor = raii::move_constructor;
      std::uint64_t old_copy_assignment = raii::copy_assignment;
      std::uint64_t old_move_assignment = raii::move_assignment;

      visitor(values, x, reference_cont);

      BOOST_TEST_EQ(old_default_constructor, raii::default_constructor);
      BOOST_TEST_EQ(old_copy_constructor, raii::copy_constructor);
      BOOST_TEST_EQ(old_move_constructor, raii::move_constructor);
      BOOST_TEST_EQ(old_copy_assignment, raii::copy_assignment);
      BOOST_TEST_EQ(old_move_assignment, raii::move_assignment);
    }

    BOOST_TEST_GE(raii::default_constructor, 0u);
    BOOST_TEST_GE(raii::copy_constructor, 0u);
    BOOST_TEST_GE(raii::move_constructor, 0u);
    BOOST_TEST_GT(raii::destructor, 0u);

    BOOST_TEST_EQ(raii::default_constructor + raii::copy_constructor +
                    raii::move_constructor,
      raii::destructor);
  }

  template <class X, class GF>
  void empty_visit(X*, GF gen_factory, test::random_generator rg)
  {
    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    using values_type = decltype(values);
    using span_value_type = typename values_type::value_type;

    raii::reset_counts();

    {
      X x;

      std::uint64_t old_default_constructor = raii::default_constructor;
      std::uint64_t old_copy_constructor = raii::copy_constructor;
      std::uint64_t old_move_constructor = raii::move_constructor;
      std::uint64_t old_copy_assignment = raii::copy_assignment;
      std::uint64_t old_move_assignment = raii::move_assignment;

      {
        thread_runner(values, [&x](boost::span<span_value_type> s) {
          std::atomic<std::uint64_t> num_visits{0};

          x.visit_all(
            [&num_visits](typename X::value_type const&) { ++num_visits; });
          BOOST_TEST_EQ(num_visits, 0u);

          for (auto const& val : s) {
            auto count = x.visit(get_key(val),
              [&num_visits](typename X::value_type const&) { ++num_visits; });
            BOOST_TEST_EQ(count, 0u);
          }
        });
      }

      BOOST_TEST_EQ(old_default_constructor, raii::default_constructor);
      BOOST_TEST_EQ(old_copy_constructor, raii::copy_constructor);
      BOOST_TEST_EQ(old_move_constructor, raii::move_constructor);
      BOOST_TEST_EQ(old_copy_assignment, raii::copy_assignment);
      BOOST_TEST_EQ(old_move_assignment, raii::move_assignment);
    }

    BOOST_TEST_EQ(raii::default_constructor, 0u);
    BOOST_TEST_EQ(raii::copy_constructor, 0u);
    BOOST_TEST_EQ(raii::move_constructor, 0u);
    BOOST_TEST_EQ(raii::destructor, 0u);
  }

  template <class X, class GF>
  void insert_and_visit(X*, GF gen_factory, test::random_generator rg)
  {
    // here we attempt to ensure happens-before and synchronizes-with
    // the visitation thread essentially chases the insertion one
    // we double-check unreloated loads/stores to ensure that a store is visible
    // in the visitation thread

    BOOST_TEST(rg == test::sequential);

    auto gen = gen_factory.template get<X>();
    auto const values = make_random_values(1024 * 16, [&] { return gen(rg); });

    {
      raii::reset_counts();

      X x;

      std::thread t1, t2;
      boost::compat::latch l(2);
      std::vector<std::string> strs(values.size());

      t1 = std::thread([&l, &values, &x, &strs] {
        l.arrive_and_wait();
        for (std::size_t idx = 0; idx < values.size(); ++idx) {
          strs[idx] = "rawr";
          auto const& val = values[idx];
          x.insert(val);
        }
      });

      t2 = std::thread([&l, &values, &x, &strs] {
        l.arrive_and_wait();

        for (std::size_t idx = 0; idx < values.size(); ++idx) {
          std::atomic_bool b{false};
          while (!b) {
            x.cvisit(get_key(values[idx]),
              [&b, &strs, idx, &values](typename X::value_type const& v) {
                BOOST_TEST_EQ(get_value(v), get_value(values[idx]));
                BOOST_TEST_EQ(strs[idx], "rawr");
                b = true;
              });
          }
        }
      });

      t1.join();
      t2.join();
    }
    check_raii_counts();
  }

  struct regular_key_extractor
  {
    template<typename T>
    auto operator()(const T& x) const -> decltype(get_key(x))
    {
      return get_key(x);
    }
  } regular_key_extract;

  struct transp_key_extractor
  {
    template<typename T>
    auto operator()(const T& x) const -> decltype((get_key(x).x_))
    {
      return get_key(x).x_;
    }
  } transp_key_extract;

  template <class X, class KeyExtractor, class GF>
  void bulk_visit(
    X*, KeyExtractor key_extract, GF gen_factory, test::random_generator rg)
  {
    using key_type = typename X::key_type;
    using value_type = typename X::value_type;

    // concurrent_flat_set visit is always const access
    using arg_type = typename std::conditional<
      std::is_same<key_type, value_type>::value,
      value_type const,
      value_type
    >::type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(16384 * 16, [&] { return gen(rg); });

    using values_type = decltype(values);
    using span_value_type = typename values_type::value_type;

    raii::reset_counts();

    {
      X x;
      for (auto const& v: values) {
        if (get_key(v).x_ % 3 != 0) x.insert(v);
      }
      X const& cx = x;

      std::uint64_t old_default_constructor = raii::default_constructor;
      std::uint64_t old_copy_constructor = raii::copy_constructor;
      std::uint64_t old_move_constructor = raii::move_constructor;
      std::uint64_t old_copy_assignment = raii::copy_assignment;
      std::uint64_t old_move_assignment = raii::move_assignment;

      std::atomic<std::size_t> num_visits{0};

      thread_runner(values, [&x, &cx, &num_visits, key_extract]
        (boost::span<span_value_type> s) {
        auto it = boost::make_transform_iterator(s.begin(), key_extract);

        std::size_t n = s.size(), m = 0, q = 0;

        auto found = [&it, &m](value_type const& v) {
          return std::find(
            it, it + (std::ptrdiff_t)m, get_key(v)) != it + (std::ptrdiff_t)m;
        };

        while (n) {
          if (m > n) m = n;

          switch (q % 3) {
            case 0:
              x.visit(
                it, it + (std::ptrdiff_t)m,
                [&num_visits, &found](arg_type& v) {
                  if ( found(v) ) ++num_visits;
                });
              break;
            case 1:
              cx.visit(
                it, it + (std::ptrdiff_t)m,
                [&num_visits, &found](value_type const& v) {
                  if ( found(v) ) ++num_visits;
                });
              break;
            case 2:
              cx.cvisit(
                it, it + (std::ptrdiff_t)m,
                [&num_visits, &found](value_type const& v) {
                  if ( found(v) ) ++num_visits;
                });
              break;
            default:
              break;
          }
          it += (std::ptrdiff_t)m;
          n -= m;
          ++m;
          if (m > 5*X::bulk_visit_size){
            m = 0;
            ++ q;
          }
        }
      });     

      BOOST_TEST_EQ(num_visits, x.size());

      BOOST_TEST_EQ(old_default_constructor, raii::default_constructor);
      BOOST_TEST_EQ(old_copy_constructor, raii::copy_constructor);
      BOOST_TEST_EQ(old_move_constructor, raii::move_constructor);
      BOOST_TEST_EQ(old_copy_assignment, raii::copy_assignment);
      BOOST_TEST_EQ(old_move_assignment, raii::move_assignment);
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
  boost::unordered::concurrent_node_map<raii, raii>* node_map;
  boost::unordered::concurrent_node_map<raii, raii, transp_hash,
    transp_key_equal>* transp_node_map;
  boost::unordered::concurrent_flat_set<raii>* set;
  boost::unordered::concurrent_flat_set<raii, transp_hash,
    transp_key_equal>* transp_set;
  boost::unordered::concurrent_node_set<raii>* node_set;
  boost::unordered::concurrent_node_set<raii, transp_hash,
    transp_key_equal>* transp_node_set;

  struct mutable_pair
  {
    mutable_pair(int first_ = 0, int second_ = 0):
      first{first_}, second{second_} {}

    int         first = 0;
    mutable int second = 0;
  };

  struct null_mutable_pair
  {
    operator mutable_pair() const { return {0,0}; }
  };

  struct mutable_pair_hash
  {
    using is_transparent = void;

    std::size_t operator()(const mutable_pair& x) const
    {
      return boost::hash<int>{}(x.first);
    }

    std::size_t operator()(const null_mutable_pair&) const
    {
      return boost::hash<int>{}(0);
    }
  };

  struct mutable_pair_equal_to
  {
    using is_transparent = void;

    bool operator()(const mutable_pair& x, const mutable_pair& y) const
    {
      return x.first == y.first;
    }

    bool operator()(const null_mutable_pair&, const mutable_pair& y) const
    {
      return 0 == y.first;
    }

    bool operator()(const mutable_pair& x, const null_mutable_pair&) const
    {
      return x.first == 0;
    }
  };

  template<typename F>
  void exclusive_access_for(F f)
  {
    std::atomic_int      num_started{0};
    std::atomic_int      in_visit{0};
    boost::compat::latch finish{1};

    auto bound_f = [&] {
      ++num_started;
      f([&] (const mutable_pair& x) {
        ++in_visit;
        ++x.second;
        finish.wait();
        --in_visit;
        return true;
      });
    };

    std::thread t1{bound_f}, t2{bound_f};
    while(num_started != 2) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    BOOST_TEST(in_visit <= 1);
    finish.count_down();
    t1.join();
    t2.join();
  }  

  template<class X>
  void exclusive_access_set_visit(X*)
  {
    using visit_function = std::function<void (const mutable_pair&)>;
    using returning_visit_function = std::function<bool (const mutable_pair&)>;
    X x;
    x.insert({0, 0});

    exclusive_access_for([&](visit_function f) {
      x.visit({0, 0}, f);
    });
    exclusive_access_for([&](visit_function f) {
      x.visit(null_mutable_pair{}, f);
    });

    exclusive_access_for([&](visit_function f) {
      mutable_pair a[] = {{0, 0}};
      x.visit(std::begin(a), std::end(a), f);
    });

    exclusive_access_for([&](visit_function f) {
      x.visit_all(f);
    });
    exclusive_access_for([&](returning_visit_function f) {
      x.visit_while(f);
    });

  #if defined(BOOST_UNORDERED_PARALLEL_ALGORITHMS)
    exclusive_access_for([&](visit_function f) {
      x.visit_all(std::execution::par, f);
    });
    exclusive_access_for([&](returning_visit_function f) {
      x.visit_while(std::execution::par, f);
    });
  #endif

    exclusive_access_for([&](visit_function f) {
      const mutable_pair p;
      x.insert_or_visit(p, f);
    });
    exclusive_access_for([&](visit_function f) {
      x.insert_or_visit({0,0}, f);
    });
    exclusive_access_for([&](visit_function f) {
      x.insert_or_visit(null_mutable_pair{}, f);
    });

    exclusive_access_for([&](visit_function f) {
      mutable_pair a[] = {{0, 0}};
      x.insert_or_visit(std::begin(a), std::end(a), f);
    });
    exclusive_access_for([&](visit_function f) {
      std::initializer_list<mutable_pair> il = {{0, 0}};
      x.insert_or_visit(il, f);
    });

    exclusive_access_for([&](visit_function f) {
      x.emplace_or_visit(0, 0, f);
    });
  }

  boost::concurrent_flat_set<
    mutable_pair, mutable_pair_hash, mutable_pair_equal_to>* mutable_set;
  boost::concurrent_node_set<
    mutable_pair, mutable_pair_hash, mutable_pair_equal_to>* mutable_node_set;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off

UNORDERED_TEST(
  visit,
  ((map)(node_map)(set)(node_set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((lvalue_visitor)(visit_all)(visit_while)(exec_policy_visit_all)
   (exec_policy_visit_while))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  visit,
  ((transp_map)(transp_node_map)(transp_set)(transp_node_set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((transp_visitor))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  empty_visit,
  ((map)(transp_map)(node_map)(transp_node_map)
   (set)(transp_set)(node_set)(transp_node_set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((default_generator)(sequential)(limited_range))
)

UNORDERED_TEST(
  insert_and_visit,
  ((map)(node_map)(set)(node_set))
  ((value_type_generator_factory))
  ((sequential))
)

UNORDERED_TEST(
  bulk_visit,
  ((map)(node_map)(set)(node_set))
  ((regular_key_extract))
  ((value_type_generator_factory))
  ((sequential))
)

UNORDERED_TEST(
  bulk_visit,
  ((transp_map)(transp_node_map)(transp_set)(transp_node_set))
  ((transp_key_extract))
  ((value_type_generator_factory))
  ((sequential))
)

// https://github.com/boostorg/unordered/issues/260

UNORDERED_TEST(
  exclusive_access_set_visit,
  ((mutable_set)(mutable_node_set))
)

// clang-format on

RUN_TESTS()
