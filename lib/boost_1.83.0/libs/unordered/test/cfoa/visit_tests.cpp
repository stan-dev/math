// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

#include <boost/core/ignore_unused.hpp>

#include <functional>
#include <vector>

namespace {
  test::seed_t initialize_seed(335740237);

  struct lvalue_visitor_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_map)
    {
      using value_type = typename X::value_type;

      std::atomic<std::uint64_t> num_visits{0};
      std::atomic<std::uint64_t> total_count{0};

      auto mut_visitor = [&num_visits, &reference_map](value_type& v) {
        BOOST_TEST(reference_map.contains(v.first));
        BOOST_TEST_EQ(v.second, reference_map.find(v.first)->second);
        ++num_visits;
      };

      auto const_visitor = [&num_visits, &reference_map](value_type const& v) {
        BOOST_TEST(reference_map.contains(v.first));
        BOOST_TEST_EQ(v.second, reference_map.find(v.first)->second);
        ++num_visits;
      };

      {
        thread_runner(
          values, [&x, &mut_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = val.first.x_;
              BOOST_TEST(r >= 0);

              auto count = x.visit(val.first, mut_visitor);
              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.visit(val.second, mut_visitor);
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
              auto r = val.first.x_;
              BOOST_TEST(r >= 0);

              auto const& y = x;
              auto count = y.visit(val.first, const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = y.visit(val.second, const_visitor);
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
              auto r = val.first.x_;
              BOOST_TEST(r >= 0);

              auto count = x.cvisit(val.first, const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.cvisit(val.second, const_visitor);
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
            auto r = val.first.x_;
            BOOST_TEST(r >= 0);

            auto count = x.count(val.first);
            BOOST_TEST_EQ(count, 1u);
            total_count += count;

            count = x.count(val.second);
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
            auto r = val.first.x_;
            BOOST_TEST(r >= 0);

            auto contains = x.contains(val.first);
            BOOST_TEST(contains);

            contains = x.contains(val.second);
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
    void operator()(std::vector<T>& values, X& x, M const& reference_map)
    {
      using value_type = typename X::value_type;

      std::atomic<std::uint64_t> num_visits{0};
      std::atomic<std::uint64_t> total_count{0};

      auto mut_visitor = [&num_visits, &reference_map](value_type& v) {
        BOOST_TEST(reference_map.contains(v.first));
        BOOST_TEST_EQ(v.second, reference_map.find(v.first)->second);
        ++num_visits;
      };

      auto const_visitor = [&num_visits, &reference_map](value_type const& v) {
        BOOST_TEST(reference_map.contains(v.first));
        BOOST_TEST_EQ(v.second, reference_map.find(v.first)->second);
        ++num_visits;
      };

      {
        thread_runner(
          values, [&x, &mut_visitor, &total_count](boost::span<T> s) {
            for (auto const& val : s) {
              auto r = val.first.x_;
              BOOST_TEST(r >= 0);

              auto count = x.visit(val.first.x_, mut_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.visit(val.second.x_, mut_visitor);
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
              auto r = val.first.x_;
              BOOST_TEST(r >= 0);

              auto const& y = x;
              auto count = y.visit(val.first.x_, const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = y.visit(val.second.x_, const_visitor);
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
              auto r = val.first.x_;
              BOOST_TEST(r >= 0);

              auto count = x.cvisit(val.first.x_, const_visitor);

              BOOST_TEST_EQ(count, 1u);
              total_count += count;

              count = x.cvisit(val.second.x_, const_visitor);
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
            auto r = val.first.x_;
            BOOST_TEST(r >= 0);

            auto count = x.count(val.first.x_);
            BOOST_TEST_EQ(count, 1u);
            total_count += count;

            count = x.count(val.second.x_);
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
            auto r = val.first.x_;
            BOOST_TEST(r >= 0);

            auto contains = x.contains(val.first.x_);
            BOOST_TEST(contains);

            contains = x.contains(val.second.x_);
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
    void operator()(std::vector<T>& values, X& x, M const& reference_map)
    {
      using value_type = typename X::value_type;

      std::atomic<std::uint64_t> total_count{0};

      auto mut_visitor = [&reference_map](std::atomic<uint64_t>& num_visits) {
        return [&reference_map, &num_visits](value_type& kv) {
          BOOST_TEST(reference_map.contains(kv.first));
          BOOST_TEST_EQ(kv.second, reference_map.find(kv.first)->second);
          ++num_visits;
        };
      };

      auto const_visitor = [&reference_map](std::atomic<uint64_t>& num_visits) {
        return [&reference_map, &num_visits](value_type const& kv) {
          BOOST_TEST(reference_map.contains(kv.first));
          BOOST_TEST_EQ(kv.second, reference_map.find(kv.first)->second);
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

  struct exec_policy_visit_all_type
  {
    template <class T, class X, class M>
    void operator()(std::vector<T>& values, X& x, M const& reference_map)
    {
#if defined(BOOST_UNORDERED_PARALLEL_ALGORITHMS)
      using value_type = typename X::value_type;

      auto mut_visitor = [&reference_map](std::atomic<uint64_t>& num_visits) {
        return [&reference_map, &num_visits](value_type& kv) {
          BOOST_TEST(reference_map.contains(kv.first));
          BOOST_TEST_EQ(kv.second, reference_map.find(kv.first)->second);
          ++num_visits;
        };
      };

      auto const_visitor = [&reference_map](std::atomic<uint64_t>& num_visits) {
        return [&reference_map, &num_visits](value_type const& kv) {
          BOOST_TEST(reference_map.contains(kv.first));
          BOOST_TEST_EQ(kv.second, reference_map.find(kv.first)->second);
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
      (void)reference_map;
#endif
    }
  } exec_policy_visit_all;

  template <class X, class G, class F>
  void visit(X*, G gen, F visitor, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    for (auto& val : values) {
      if (val.second.x_ == 0) {
        val.second.x_ = 1;
      }
      val.second.x_ *= -1;
    }

    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());

    raii::reset_counts();

    {
      X x;
      for (auto const& v : values) {
        x.insert(v);
      }
      BOOST_TEST_EQ(x.size(), reference_map.size());

      std::uint64_t old_default_constructor = raii::default_constructor;
      std::uint64_t old_copy_constructor = raii::copy_constructor;
      std::uint64_t old_move_constructor = raii::move_constructor;
      std::uint64_t old_copy_assignment = raii::copy_assignment;
      std::uint64_t old_move_assignment = raii::move_assignment;

      visitor(values, x, reference_map);

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

  template <class X, class G>
  void empty_visit(X*, G gen, test::random_generator rg)
  {
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
            auto count = x.visit(val.first,
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

  template <class X, class G>
  void insert_and_visit(X*, G gen, test::random_generator rg)
  {
    // here we attempt to ensure happens-before and synchronizes-with
    // the visitation thread essentially chases the insertion one
    // we double-check unreloated loads/stores to ensure that a store is visible
    // in the visitation thread

    BOOST_TEST(rg == test::sequential);

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
            x.cvisit(values[idx].first,
              [&b, &strs, idx, &values](typename X::value_type const& v) {
                BOOST_TEST_EQ(v.second, values[idx].second);
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

  boost::unordered::concurrent_flat_map<raii, raii>* map;
  boost::unordered::concurrent_flat_map<raii, raii, transp_hash,
    transp_key_equal>* transp_map;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off

UNORDERED_TEST(
  visit,
  ((map))
  ((value_type_generator)(init_type_generator))
  ((lvalue_visitor)(visit_all)(exec_policy_visit_all))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  visit,
  ((transp_map))
  ((value_type_generator)(init_type_generator))
  ((transp_visitor))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  empty_visit,
  ((map)(transp_map))
  ((value_type_generator)(init_type_generator))
  ((default_generator)(sequential)(limited_range))
)

UNORDERED_TEST(
  insert_and_visit,
  ((map))
  ((value_type_generator))
  ((sequential))
)

// clang-format on

RUN_TESTS()
