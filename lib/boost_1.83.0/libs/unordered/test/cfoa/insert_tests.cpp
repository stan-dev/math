// Copyright (C) 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>

#include <boost/core/ignore_unused.hpp>

struct raii_convertible
{
  int x, y;
  raii_convertible(int x_, int y_) : x{x_}, y{y_} {}

  operator std::pair<raii const, raii>() { return {x, y}; }
};

namespace {
  test::seed_t initialize_seed(78937);

  struct lvalue_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
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
      BOOST_TEST_EQ(raii::copy_constructor, 2 * x.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_inserter;

  struct norehash_lvalue_inserter_type : public lvalue_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      x.reserve(values.size());
      lvalue_inserter_type::operator()(values, x);
      BOOST_TEST_EQ(raii::copy_constructor, 2 * x.size());
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

      if (std::is_same<T, typename X::value_type>::value) {
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
      x.reserve(values.size());

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_EQ(raii::move_constructor, 0u);

      rvalue_inserter_type::operator()(values, x);

      if (std::is_same<T, typename X::value_type>::value) {
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
        BOOST_TEST_EQ(raii::move_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_EQ(raii::move_constructor, 2 * x.size());
      }
    }
  } norehash_rvalue_inserter;

  struct iterator_range_inserter_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& p : values) {
        values2.push_back(raii_convertible(p.first.x_, p.second.x_));
      }

      thread_runner(values2, [&x](boost::span<raii_convertible> s) {
        x.insert(s.begin(), s.end());
      });

      BOOST_TEST_EQ(raii::default_constructor, 2 * values2.size());
#if BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
    BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)
      // some versions of old gcc have trouble eliding copies here
      // https://godbolt.org/z/Ebo6TbvaG
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
      // don't check move construction count here because of rehashing
      BOOST_TEST_GT(raii::move_constructor, 0u);
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
      BOOST_TEST_GT(raii::move_constructor, x.size()); // rehashing
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
      BOOST_TEST_GT(raii::move_constructor, x.size()); // rehashing
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
      BOOST_TEST_GE(raii::move_constructor, 2 * x.size());
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
      BOOST_TEST_GT(raii::move_constructor, x.size()); // rehashing
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
      BOOST_TEST_GT(raii::move_constructor, 2 * x.size()); // rehashing
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, values.size() - x.size());
    }
  } trans_insert_or_assign_move_assign;

  struct lvalue_insert_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
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
      BOOST_TEST_EQ(raii::copy_constructor, 2 * x.size());
      // don't check move construction count here because of rehashing
      BOOST_TEST_GT(raii::move_constructor, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_insert_or_cvisit;

  struct lvalue_insert_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b =
            x.insert_or_visit(r, [&num_invokes](typename X::value_type& v) {
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
      BOOST_TEST_EQ(raii::copy_constructor, 2 * x.size());
      // don't check move construction count here because of rehashing
      BOOST_TEST_GT(raii::move_constructor, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } lvalue_insert_or_visit;

  struct rvalue_insert_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
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
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
        BOOST_TEST_GE(raii::move_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(raii::move_constructor, 2 * x.size());
      }
    }
  } rvalue_insert_or_cvisit;

  struct rvalue_insert_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::atomic<std::uint64_t> num_inserts{0};
      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(values, [&x, &num_inserts, &num_invokes](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.insert_or_visit(
            std::move(r), [&num_invokes](typename X::value_type& v) {
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
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
        BOOST_TEST_GE(raii::move_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
        BOOST_TEST_GE(raii::move_constructor, 2 * x.size());
      }
    }
  } rvalue_insert_or_visit;

  struct iterator_range_insert_or_cvisit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& p : values) {
        values2.push_back(raii_convertible(p.first.x_, p.second.x_));
      }

      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(
        values2, [&x, &num_invokes](boost::span<raii_convertible> s) {
          x.insert_or_cvisit(s.begin(), s.end(),
            [&num_invokes](typename X::value_type const& v) {
              (void)v;
              ++num_invokes;
            });
        });

      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 2 * values2.size());
#if BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
    BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)
      // skip test
#else
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
#endif
      BOOST_TEST_GT(raii::move_constructor, 0u);
    }
  } iterator_range_insert_or_cvisit;

  struct iterator_range_insert_or_visit_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      std::vector<raii_convertible> values2;
      values2.reserve(values.size());
      for (auto const& p : values) {
        values2.push_back(raii_convertible(p.first.x_, p.second.x_));
      }

      std::atomic<std::uint64_t> num_invokes{0};
      thread_runner(
        values2, [&x, &num_invokes](boost::span<raii_convertible> s) {
          x.insert_or_visit(s.begin(), s.end(),
            [&num_invokes](typename X::value_type const& v) {
              (void)v;
              ++num_invokes;
            });
        });

      BOOST_TEST_EQ(num_invokes, values.size() - x.size());

      BOOST_TEST_EQ(raii::default_constructor, 2 * values2.size());
#if BOOST_WORKAROUND(BOOST_GCC_VERSION, >= 50300) && \
    BOOST_WORKAROUND(BOOST_GCC_VERSION, <  50500)
      // skip test
#else
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
#endif
      BOOST_TEST_GT(raii::move_constructor, 0u);
    }
  } iterator_range_insert_or_visit;

  template <class X, class G, class F>
  void insert(X*, G gen, F inserter, test::random_generator rg)
  {
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());
    raii::reset_counts();

    {
      X x;

      inserter(values, x);

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

  template <class X> void insert_initializer_list(X*)
  {
    using value_type = typename X::value_type;

    std::initializer_list<value_type> values{
      value_type{raii{0}, raii{0}},
      value_type{raii{1}, raii{1}},
      value_type{raii{2}, raii{2}},
      value_type{raii{3}, raii{3}},
      value_type{raii{4}, raii{4}},
      value_type{raii{5}, raii{5}},
      value_type{raii{6}, raii{6}},
      value_type{raii{6}, raii{6}},
      value_type{raii{7}, raii{7}},
      value_type{raii{8}, raii{8}},
      value_type{raii{9}, raii{9}},
      value_type{raii{10}, raii{10}},
      value_type{raii{9}, raii{9}},
      value_type{raii{8}, raii{8}},
      value_type{raii{7}, raii{7}},
      value_type{raii{6}, raii{6}},
      value_type{raii{5}, raii{5}},
      value_type{raii{4}, raii{4}},
      value_type{raii{3}, raii{3}},
      value_type{raii{2}, raii{2}},
      value_type{raii{1}, raii{1}},
      value_type{raii{0}, raii{0}},
    };

    std::vector<raii> dummy;

    auto reference_map =
      boost::unordered_flat_map<raii, raii>(values.begin(), values.end());
    raii::reset_counts();

    {
      {
        X x;

        thread_runner(
          dummy, [&x, &values](boost::span<raii>) { x.insert(values); });

        BOOST_TEST_EQ(x.size(), reference_map.size());

        BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& kv) {
          BOOST_TEST(reference_map.contains(kv.first));
          BOOST_TEST_EQ(kv.second, reference_map[kv.first]);
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

        thread_runner(dummy, [&x, &values, &num_invokes](boost::span<raii>) {
          x.insert_or_visit(values, [&num_invokes](typename X::value_type& v) {
            (void)v;
            ++num_invokes;
          });

          x.insert_or_cvisit(
            values, [&num_invokes](typename X::value_type const& v) {
              (void)v;
              ++num_invokes;
            });
        });

        BOOST_TEST_EQ(num_invokes, (values.size() - x.size()) +
                                     (num_threads - 1) * values.size() +
                                     num_threads * values.size());
        BOOST_TEST_EQ(x.size(), reference_map.size());

        BOOST_TEST_EQ(x.size(), x.visit_all([&](value_type const& kv) {
          BOOST_TEST(reference_map.contains(kv.first));
          BOOST_TEST_EQ(kv.second, reference_map[kv.first]);
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

  UNORDERED_AUTO_TEST (insert_sfinae_test) {
    // mostly a compile-time tests to ensure that there's no ambiguity when a
    // user does this
    using value_type =
      typename boost::unordered::concurrent_flat_map<raii, raii>::value_type;
    boost::unordered::concurrent_flat_map<raii, raii> x;
    x.insert({1, 2});

    x.insert_or_visit({2, 3}, [](value_type&) {});
    x.insert_or_cvisit({3, 4}, [](value_type const&) {});
  }

  boost::unordered::concurrent_flat_map<raii, raii>* map;
  boost::unordered::concurrent_flat_map<raii, raii, transp_hash,
    transp_key_equal>* trans_map;
  boost::unordered::concurrent_flat_map<raii, raii, boost::hash<raii>,
    std::equal_to<raii>, fancy_allocator<std::pair<raii const, raii> > >*
    fancy_map;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off
UNORDERED_TEST(
  insert_initializer_list,
  ((map)))

UNORDERED_TEST(
  insert,
  ((map)(fancy_map))
  ((value_type_generator)(init_type_generator))
  ((lvalue_inserter)(rvalue_inserter)(iterator_range_inserter)
   (norehash_lvalue_inserter)(norehash_rvalue_inserter)
   (lvalue_insert_or_cvisit)(lvalue_insert_or_visit)
   (rvalue_insert_or_cvisit)(rvalue_insert_or_visit)
   (iterator_range_insert_or_cvisit)(iterator_range_insert_or_visit))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  insert,
  ((map))
  ((init_type_generator))
  ((lvalue_insert_or_assign_copy_assign)(lvalue_insert_or_assign_move_assign)
   (rvalue_insert_or_assign_copy_assign)(rvalue_insert_or_assign_move_assign))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  insert,
  ((trans_map))
  ((init_type_generator))
  ((trans_insert_or_assign_copy_assign)(trans_insert_or_assign_move_assign))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
