// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Copyright (C) 2024 Braden Ganetsky
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"
#include "../helpers/count.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

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

  template <typename Container, typename Value, typename F1, typename F2>
  bool member_emplace_and_visit(Container& x, Value& v, F1 f1, F2 f2)
  {
    return x.emplace_and_visit(v.x_, f1, f2);
  }

  template <
    typename Container, typename Key, typename Value, typename F1, typename F2>
  bool member_emplace_and_visit(
    Container& x, std::pair<Key, Value>& v, F1 f1, F2 f2)
  {
    return x.emplace_and_visit(v.first.x_, v.second.x_, f1, f2);
  }

  template <typename Container, typename Value, typename F1, typename F2>
  bool member_emplace_and_cvisit(Container& x, Value& v, F1 f1, F2 f2)
  {
    return x.emplace_and_cvisit(v.x_, f1, f2);
  }

  template <
    typename Container, typename Key, typename Value, typename F1, typename F2>
  bool member_emplace_and_cvisit(
    Container& x, std::pair<Key, Value>& v, F1 f1, F2 f2)
  {
    return x.emplace_and_cvisit(v.first.x_, v.second.x_, f1, f2);
  }

  struct lvalue_emplacer_type
  {
    template <class T, class X> void call_impl(std::vector<T>& values, X& x)
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

      std::uint64_t const default_constructors = value_type_cardinality == 2
                                                   ? values.size() + num_inserts
                                                   : values.size();
      BOOST_TEST_EQ(raii::default_constructor, default_constructors);

      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      call_impl(values, x);
      BOOST_TEST_GE(raii::move_constructor, x.size());
    }
  } lvalue_emplacer;

  struct norehash_lvalue_emplacer_type : public lvalue_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      x.reserve(values.size());
      lvalue_emplacer_type::call_impl(values, x);
      BOOST_TEST_EQ(raii::move_constructor, x.size());
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

  struct lvalue_emplace_and_cvisit_type
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
          bool b = member_emplace_and_cvisit(
            x, r,
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

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_emplace_and_cvisit;

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

  struct lvalue_emplace_and_visit_type
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
          bool b = member_emplace_and_visit(
            x, r,
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

      BOOST_TEST_EQ(
        raii::default_constructor, value_type_cardinality * values.size());
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
      BOOST_TEST_GE(raii::move_constructor, value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
    }
  } lvalue_emplace_and_visit;

  struct copy_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto value_type_cardinality =
        value_cardinality<typename X::value_type>::value;

      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto const& r : s) {
          bool b = x.emplace(r);
          if (b) {
            ++num_inserts;
          }
        }
      });
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(raii::default_constructor, 0u);

      BOOST_TEST_EQ(raii::copy_constructor, value_type_cardinality * x.size());
      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(raii::move_constructor, 0u);
      }
      else {
        BOOST_TEST_GT(raii::move_constructor, 0u);
      }

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } copy_emplacer;

  struct move_emplacer_type
  {
    template <class T, class X> void operator()(std::vector<T>& values, X& x)
    {
      static constexpr auto input_type_nonconst_cardinality = 
        value_nonconst_cardinality<T>::value;

      std::atomic<std::uint64_t> num_inserts{0};
      thread_runner(values, [&x, &num_inserts](boost::span<T> s) {
        for (auto& r : s) {
          bool b = x.emplace(std::move(r));
          if (b) {
            ++num_inserts;
          }
        }
      });
      BOOST_TEST_EQ(num_inserts, x.size());
      BOOST_TEST_EQ(raii::default_constructor, 0u);

#if defined(BOOST_MSVC)
#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant
#endif
      if (std::is_same<T, typename X::value_type>::value &&
          !std::is_same<typename X::key_type,
            typename X::value_type>::value) { // map value_type can only be
                                              // copied, no move
        BOOST_TEST_EQ(raii::copy_constructor, x.size());
      } else {
        BOOST_TEST_EQ(raii::copy_constructor, 0u);
      }

      if (is_container_node_based<X>::value) {
        BOOST_TEST_EQ(
          raii::move_constructor, input_type_nonconst_cardinality * x.size());
      }
      else {
        BOOST_TEST_GT(
          raii::move_constructor, input_type_nonconst_cardinality * x.size());
      }
#if defined(BOOST_MSVC)
#pragma warning(pop) // C4127
#endif

      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
  } move_emplacer;

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
  boost::unordered::concurrent_node_map<raii, raii>* node_map;
  boost::unordered::concurrent_flat_set<raii>* set;
  boost::unordered::concurrent_node_set<raii>* node_set;

} // namespace

using test::default_generator;
using test::limited_range;
using test::sequential;

// clang-format off

UNORDERED_TEST(
  emplace,
  ((map)(node_map)(set)(node_set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((lvalue_emplacer)(norehash_lvalue_emplacer)
   (lvalue_emplace_or_cvisit)(lvalue_emplace_or_visit)(copy_emplacer)(move_emplacer))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  emplace,
  ((map)(node_map)(set)(node_set))
  ((value_type_generator_factory)(init_type_generator_factory))
  ((lvalue_emplace_and_cvisit)(lvalue_emplace_and_visit))
  ((default_generator)(sequential)(limited_range)))

// clang-format on

namespace {
  using converting_key_type = basic_raii<struct converting_key_tag_>;
  using converting_value_type = basic_raii<struct converting_value_tag_>;

  class counted_key_type : public basic_raii<struct counted_key_tag_>
  {
  public:
    using basic_raii::basic_raii;
    counted_key_type() = default;
    counted_key_type(const converting_key_type& k) : counted_key_type(k.x_) {}
  };
  class counted_value_type : public basic_raii<struct counted_value_tag_>
  {
  public:
    using basic_raii::basic_raii;
    counted_value_type() = default;
    counted_value_type(const converting_value_type& v)
        : counted_value_type(v.x_)
    {
    }
  };

  void reset_counts()
  {
    counted_key_type::reset_counts();
    counted_value_type::reset_counts();
    converting_key_type::reset_counts();
    converting_value_type::reset_counts();
  }

  using test::smf_count;

  template <class T> smf_count count_for()
  {
    return test::smf_count{
      (int)T::default_constructor.load(std::memory_order_relaxed),
      (int)T::copy_constructor.load(std::memory_order_relaxed),
      (int)T::move_constructor.load(std::memory_order_relaxed),
      (int)T::copy_assignment.load(std::memory_order_relaxed),
      (int)T::move_assignment.load(std::memory_order_relaxed),
      (int)T::destructor.load(std::memory_order_relaxed)};
  }

  enum emplace_kind
  {
    copy,
    move
  };

  enum emplace_status
  {
    fail,
    success
  };

  struct counted_key_checker_type
  {
    using key_type = counted_key_type;
    void operator()(emplace_kind kind, emplace_status status)
    {
      int copies = (kind == copy && status == success) ? 1 : 0;
      int moves = (kind == move && status == success) ? 1 : 0;
      BOOST_TEST_EQ(
        count_for<counted_key_type>(), (smf_count{0, copies, moves, 0, 0, 0}));
    }
  } counted_key_checker;

  struct converting_key_checker_type
  {
    using key_type = converting_key_type;
    void operator()(emplace_kind, emplace_status status)
    {
      int moves = (status == success) ? 1 : 0;
      BOOST_TEST_EQ(
        count_for<counted_key_type>(), (smf_count{1, 0, moves, 0, 0, 1}));
    }
  } converting_key_checker;

  struct counted_value_checker_type
  {
    using mapped_type = counted_value_type;
    void operator()(emplace_kind kind, emplace_status status)
    {
      int copies = (kind == copy && status == success) ? 1 : 0;
      int moves = (kind == move && status == success) ? 1 : 0;
      BOOST_TEST_EQ(count_for<counted_value_type>(),
        (smf_count{0, copies, moves, 0, 0, 0}));
    }
  } counted_value_checker;

  struct converting_value_checker_type
  {
    using mapped_type = converting_value_type;
    void operator()(emplace_kind, emplace_status status)
    {
      int ctors = (status == success) ? 1 : 0;
      BOOST_TEST_EQ(
        count_for<counted_value_type>(), (smf_count{ctors, 0, 0, 0, 0, 0}));
    }
  } converting_value_checker;

  template <class X, class KC, class VC>
  void emplace_map_key_value(
    X*, emplace_kind kind, KC key_checker, VC value_checker)
  {
    using container = X;
    using key_type = typename KC::key_type;
    using mapped_type = typename VC::mapped_type;

    container x;
    key_type key{};
    key_type key2 = key;
    mapped_type value{};
    mapped_type value2 = value;

    {
      reset_counts();
      auto ret = (kind == copy) ? x.emplace(key, value)
                                : x.emplace(std::move(key), std::move(value));
      BOOST_TEST_EQ(ret, true);
      key_checker(kind, success);
      value_checker(kind, success);
      BOOST_TEST_EQ(
        count_for<converting_key_type>(), (smf_count{0, 0, 0, 0, 0, 0}));
      BOOST_TEST_EQ(
        count_for<converting_value_type>(), (smf_count{0, 0, 0, 0, 0, 0}));
    }

    {
      reset_counts();
      bool ret = x.emplace(key2, value2);
      BOOST_TEST_EQ(ret, false);
      key_checker(kind, fail);
      value_checker(kind, fail);
      BOOST_TEST_EQ(
        count_for<converting_key_type>(), (smf_count{0, 0, 0, 0, 0, 0}));
      BOOST_TEST_EQ(
        count_for<converting_value_type>(), (smf_count{0, 0, 0, 0, 0, 0}));
    }

    {
      reset_counts();
      bool ret = x.emplace(std::move(key2), std::move(value2));
      BOOST_TEST_EQ(ret, false);
      key_checker(kind, fail);
      value_checker(kind, fail);
      BOOST_TEST_EQ(
        count_for<converting_key_type>(), (smf_count{0, 0, 0, 0, 0, 0}));
      BOOST_TEST_EQ(
        count_for<converting_value_type>(), (smf_count{0, 0, 0, 0, 0, 0}));
    }
  }

  boost::unordered::concurrent_flat_map<counted_key_type, counted_value_type>*
    test_counted_flat_map = {};
  boost::unordered::concurrent_node_map<counted_key_type, counted_value_type>*
    test_counted_node_map = {};

} // namespace

// clang-format off

UNORDERED_TEST(
  emplace_map_key_value,
  ((test_counted_flat_map)(test_counted_node_map))
  ((copy)(move))
  ((counted_key_checker)(converting_key_checker))
  ((counted_value_checker)(converting_value_checker))
)

// clang-format on

RUN_TESTS()
