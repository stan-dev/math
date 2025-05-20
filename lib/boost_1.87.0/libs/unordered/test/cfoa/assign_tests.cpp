// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include "../helpers/replace_allocator.hpp"
#include "../objects/non_default_ctble_allocator.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

#include <vector>

#if defined(__clang__) && defined(__has_warning)

#if __has_warning("-Wself-assign-overloaded")
#pragma clang diagnostic ignored "-Wself-assign-overloaded"
#endif

#if __has_warning("-Wself-move")
#pragma clang diagnostic ignored "-Wself-move"
#endif

#endif /* defined(__clang__) && defined(__has_warning) */

#if defined(BOOST_GCC) && BOOST_GCC >= 130000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-move"
#endif

test::seed_t initialize_seed{2762556623};

using test::default_generator;
using test::limited_range;
using test::sequential;

using hasher = stateful_hash;
using key_equal = stateful_key_equal;

using map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, stateful_allocator<std::pair<raii const, raii> > >;

using node_map_type = boost::unordered::concurrent_node_map<raii, raii, hasher,
  key_equal, stateful_allocator<std::pair<raii const, raii> > >;

using set_type = boost::unordered::concurrent_flat_set<raii, hasher,
  key_equal, stateful_allocator<raii> >;

using node_set_type = boost::unordered::concurrent_node_set<raii, hasher,
  key_equal, stateful_allocator<raii> >;

using fancy_map_type = boost::unordered::concurrent_flat_map<raii, raii, hasher,
  key_equal, stateful_allocator2<std::pair<raii const, raii> > >;

using fancy_node_map_type = boost::unordered::concurrent_node_map<raii, raii, hasher,
  key_equal, stateful_allocator2<std::pair<raii const, raii> > >;

using fancy_set_type = boost::unordered::concurrent_flat_set<raii, hasher,
  key_equal, stateful_allocator2<raii> >;

using fancy_node_set_type = boost::unordered::concurrent_node_set<raii, hasher,
  key_equal, stateful_allocator2<raii> >;

map_type* test_map;
node_map_type* test_node_map;
set_type* test_set;
node_set_type* test_node_set;
fancy_map_type* fancy_test_map;
fancy_node_map_type* fancy_test_node_map;
fancy_set_type* fancy_test_set;
fancy_node_set_type* fancy_test_node_set;

std::initializer_list<map_type::value_type> map_init_list{
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

std::initializer_list<set_type::value_type> set_init_list{
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

auto test_map_and_init_list=std::make_pair(test_map,map_init_list);
auto test_node_map_and_init_list=std::make_pair(test_node_map,map_init_list);
auto test_set_and_init_list=std::make_pair(test_set,set_init_list);
auto test_node_set_and_init_list=std::make_pair(test_node_set,set_init_list);

template <class T,bool POCCA, bool POCMA>
struct poca_allocator: fancy_allocator<T>
{
  using super = fancy_allocator<T>;
  using pointer = typename super::pointer;
  using propagate_on_container_copy_assignment = 
    std::integral_constant<bool, POCCA>;
  using propagate_on_container_move_assignment = 
    std::integral_constant<bool, POCMA>;

  int x_ = -1;

  template <class U> struct rebind
  {
    typedef poca_allocator<U, POCCA, POCMA> other;
  };

  poca_allocator() = default;
  poca_allocator(poca_allocator const&) = default;
  poca_allocator(poca_allocator &&) = default;

  poca_allocator(int const x) : x_{x} {}

  poca_allocator& operator=(poca_allocator const& rhs)
  {
    if (this != &rhs) {
      super::operator=(rhs);
      x_ = rhs.x_;
    }
    return *this;
  }

  template <class U> poca_allocator(
    poca_allocator<U, POCCA, POCMA> const& rhs) :
    super{rhs}, x_{rhs.x_}
  {
  }

  pointer allocate(std::size_t n)
  {
    auto p = super::allocate(n + 1);
    reinterpret_cast<char&>(*p) = static_cast<char>(x_);
    return p + std::ptrdiff_t(1);
  }

  void deallocate(pointer p, std::size_t n)
  {
    p = p + std::ptrdiff_t(-1);
    BOOST_TEST_EQ(reinterpret_cast<char&>(*p), static_cast<char>(x_));
    super::deallocate(p, n + 1); 
  }

  bool operator==(poca_allocator const& rhs) const { return x_ == rhs.x_; }
  bool operator!=(poca_allocator const& rhs) const { return x_ != rhs.x_; }
};

template <class T> 
struct pocca_allocator: poca_allocator<T, true, false>
{
  pocca_allocator() = default;
  pocca_allocator(pocca_allocator const&) = default;
  pocca_allocator(pocca_allocator &&) = default;
  using poca_allocator<T, true, false>::poca_allocator;

  pocca_allocator& operator=(pocca_allocator const&) = default;
};

template <class T> 
struct pocma_allocator: poca_allocator<T, false, true>
{
  pocma_allocator() = default;
  pocma_allocator(pocma_allocator const&) = default;
  pocma_allocator(pocma_allocator &&) = default;
  using poca_allocator<T, false, true>::poca_allocator;

  pocma_allocator& operator=(pocma_allocator const&) = default;
};

namespace {
  template <class X, class GF>
  void copy_assign(X*, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    // lhs empty, rhs empty
    {
      raii::reset_counts();

      X x(0, hasher(1), key_equal(2), allocator_type(3));

      thread_runner(values, [&x](boost::span<value_type> s) {
        (void)s;

        X y;

        BOOST_TEST(x.empty());
        BOOST_TEST(y.empty());

        y = x;

        BOOST_TEST_EQ(x.hash_function(), y.hash_function());
        BOOST_TEST_EQ(x.key_eq(), y.key_eq());
        BOOST_TEST(x.get_allocator() != y.get_allocator());
      });

      BOOST_TEST_EQ(raii::destructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
    }

    // lhs non-empty, rhs empty
    {
      raii::reset_counts();

      X x(0, hasher(1), key_equal(2), allocator_type(3));

      auto const old_size = reference_cont.size();

      thread_runner(values, [&x, &values](boost::span<value_type> s) {
        (void)s;

        X y(values.size());
        for (auto const& v : values) {
          y.insert(v);
        }

        BOOST_TEST(x.empty());
        BOOST_TEST(!y.empty());

        y = x;

        BOOST_TEST_EQ(x.hash_function(), y.hash_function());
        BOOST_TEST_EQ(x.key_eq(), y.key_eq());
        BOOST_TEST(x.get_allocator() != y.get_allocator());

        BOOST_TEST(y.empty());
      });

      BOOST_TEST_EQ(
        raii::destructor, num_threads * (value_type_cardinality * old_size));
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        num_threads * value_type_cardinality * reference_cont.size());
    }
    check_raii_counts();

    // lhs empty, rhs non-empty
    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_cc = +raii::copy_constructor;

      thread_runner(
        values, [&x, &reference_cont](boost::span<value_type> s) {
          (void)s;

          X y;

          BOOST_TEST(!x.empty());
          BOOST_TEST(y.empty());

          y = x;

          BOOST_TEST_EQ(x.hash_function(), y.hash_function());
          BOOST_TEST_EQ(x.key_eq(), y.key_eq());
          BOOST_TEST(x.get_allocator() != y.get_allocator());

          test_matches_reference(y, reference_cont);
        });

      BOOST_TEST_EQ(
        raii::destructor, num_threads * value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        old_cc + (num_threads * value_type_cardinality * x.size()));
    }
    check_raii_counts();

    // lhs non-empty, rhs non-empty
    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_size = x.size();
      auto const old_cc = +raii::copy_constructor;

      thread_runner(values, [&x, &values](boost::span<value_type> s) {
        (void)s;

        X y(values.size());
        for (auto const& v : values) {
          y.insert(v);
        }

        BOOST_TEST(!x.empty());
        BOOST_TEST(!y.empty());

        y = x;

        BOOST_TEST_EQ(x.hash_function(), y.hash_function());
        BOOST_TEST_EQ(x.key_eq(), y.key_eq());
        BOOST_TEST(x.get_allocator() != y.get_allocator());
      });

      BOOST_TEST_EQ(
        raii::destructor, 2 * num_threads * value_type_cardinality * old_size);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        old_cc + (2 * num_threads * value_type_cardinality * x.size()));
    }
    check_raii_counts();

    // self-assign
    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_cc = +raii::copy_constructor;

      thread_runner(
        values, [&x, &reference_cont](boost::span<value_type> s) {
          (void)s;

          BOOST_TEST(!x.empty());

          x = x;

          BOOST_TEST_EQ(x.hash_function(), hasher(1));
          BOOST_TEST_EQ(x.key_eq(), key_equal(2));
          BOOST_TEST(x.get_allocator() == allocator_type(3));

          test_matches_reference(x, reference_cont);
        });

      BOOST_TEST_EQ(raii::destructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
    }
    check_raii_counts();

    // propagation
    {
      using pocca_container_type = replace_allocator<X, pocca_allocator>;
      using pocca_allocator_type =
        typename pocca_container_type::allocator_type;

      raii::reset_counts();

      pocca_container_type x(
        values.size(), hasher(1), key_equal(2), pocca_allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_size = x.size();
      auto const old_cc = +raii::copy_constructor;

      thread_runner(values, [&x, &values](boost::span<value_type> s) {
        (void)s;

        pocca_container_type y(values.size());
        for (auto const& v : values) {
          y.insert(v);
        }

        BOOST_TEST(!x.empty());
        BOOST_TEST(!y.empty());

        BOOST_TEST(x.get_allocator() != y.get_allocator());

        y = x;

        BOOST_TEST_EQ(x.hash_function(), y.hash_function());
        BOOST_TEST_EQ(x.key_eq(), y.key_eq());
        BOOST_TEST(x.get_allocator() == y.get_allocator());
      });

      BOOST_TEST_EQ(
        raii::destructor, 2 * num_threads * value_type_cardinality * old_size);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        old_cc + (2 * num_threads * value_type_cardinality * x.size()));
    }
    check_raii_counts();
  }

  template <class X, class GF>
  void move_assign(X*, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    using pocma_container_type = replace_allocator<X, pocma_allocator>;
    using pocma_allocator_type = typename pocma_container_type::allocator_type;

    auto gen = gen_factory.template get<X>();

    BOOST_STATIC_ASSERT(
      std::is_nothrow_move_assignable<
        replace_allocator<X, std::allocator> >::value);

    BOOST_STATIC_ASSERT(
      !std::is_nothrow_move_assignable<
        replace_allocator<X, stateful_allocator> >::value);

    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    // move assignment has more complex requirements than copying
    // equal allocators:
    // lhs empty, rhs non-empty
    // lhs non-empty, rhs empty
    // lhs non-empty, rhs non-empty
    //
    // unequal allocators:
    // lhs non-empty, rhs non-empty
    //
    // pocma
    // self move-assign

    // lhs empty, rhs empty
    {
      raii::reset_counts();

      X x(0, hasher(1), key_equal(2), allocator_type(3));

      std::atomic<unsigned> num_transfers{0};

      thread_runner(
        values, [&x, &num_transfers](boost::span<value_type> s) {
          (void)s;

          X y(0, hasher(2), key_equal(1), allocator_type(3));

          BOOST_TEST(x.empty());
          BOOST_TEST(y.empty());
          BOOST_TEST(x.get_allocator() == y.get_allocator());

          y = std::move(x);
          if (y.hash_function() == hasher(1)) {
            ++num_transfers;
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.hash_function(), hasher(2));
            BOOST_TEST_EQ(y.key_eq(), key_equal(1));
          }

          BOOST_TEST_EQ(x.hash_function(), hasher(2));
          BOOST_TEST_EQ(x.key_eq(), key_equal(1));
          BOOST_TEST(x.get_allocator() == y.get_allocator());
        });

      BOOST_TEST_EQ(num_transfers, 1u);

      BOOST_TEST_EQ(raii::destructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, 0u);
    }

    // lhs non-empty, rhs empty
    {
      raii::reset_counts();

      X x(0, hasher(1), key_equal(2), allocator_type(3));

      std::atomic<unsigned> num_transfers{0};

      thread_runner(
        values, [&x, &values, &num_transfers](boost::span<value_type> s) {
          (void)s;

          X y(values.size(), hasher(2), key_equal(1), allocator_type(3));
          for (auto const& v : values) {
            y.insert(v);
          }

          BOOST_TEST(x.empty());
          BOOST_TEST(!y.empty());
          BOOST_TEST(x.get_allocator() == y.get_allocator());

          y = std::move(x);
          if (y.hash_function() == hasher(1)) {
            ++num_transfers;
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.hash_function(), hasher(2));
            BOOST_TEST_EQ(y.key_eq(), key_equal(1));
          }

          BOOST_TEST_EQ(x.hash_function(), hasher(2));
          BOOST_TEST_EQ(x.key_eq(), key_equal(1));
          BOOST_TEST(x.get_allocator() == y.get_allocator());

          BOOST_TEST(y.empty());
        });

      BOOST_TEST_EQ(num_transfers, 1u);

      BOOST_TEST_EQ(
        raii::destructor, num_threads * value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        num_threads * value_type_cardinality * reference_cont.size());
    }
    check_raii_counts();

    // lhs empty, rhs non-empty
    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;
      std::atomic<unsigned> num_transfers{0};

      thread_runner(values,
        [&x, &reference_cont, &num_transfers](boost::span<value_type> s) {
          (void)s;

          X y(allocator_type(3));

          BOOST_TEST(y.empty());
          BOOST_TEST(x.get_allocator() == y.get_allocator());

          y = std::move(x);
          if (!y.empty()) {
            ++num_transfers;
            test_matches_reference(y, reference_cont);

            BOOST_TEST_EQ(y.hash_function(), hasher(1));
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.hash_function(), hasher());
            BOOST_TEST_EQ(y.key_eq(), key_equal());
          }

          BOOST_TEST(x.empty());

          BOOST_TEST_EQ(x.hash_function(), hasher());
          BOOST_TEST_EQ(x.key_eq(), key_equal());
          BOOST_TEST(x.get_allocator() == y.get_allocator());
        });

      BOOST_TEST_EQ(num_transfers, 1u);

      BOOST_TEST_EQ(
        raii::destructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
    }
    check_raii_counts();

    // lhs non-empty, rhs non-empty
    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_size = x.size();
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      std::atomic<unsigned> num_transfers{0};

      thread_runner(values, [&x, &values, &num_transfers, &reference_cont](
                              boost::span<value_type> s) {
        (void)s;

        X y(values.size(), hasher(2), key_equal(1), allocator_type(3));
        for (auto const& v : values) {
          y.insert(v);
        }

        BOOST_TEST(!y.empty());
        BOOST_TEST(x.get_allocator() == y.get_allocator());

        y = std::move(x);
        if (y.hash_function() == hasher(1)) {
          ++num_transfers;
          test_matches_reference(y, reference_cont);

          BOOST_TEST_EQ(y.key_eq(), key_equal(2));
        } else {
          BOOST_TEST_EQ(y.hash_function(), hasher(2));
          BOOST_TEST_EQ(y.key_eq(), key_equal(1));
        }

        BOOST_TEST(x.empty());

        BOOST_TEST_EQ(x.hash_function(), hasher(2));
        BOOST_TEST_EQ(x.key_eq(), key_equal(1));
        BOOST_TEST(x.get_allocator() == y.get_allocator());
      });

      BOOST_TEST_EQ(num_transfers, 1u);

      BOOST_TEST_EQ(
        raii::destructor,
        value_type_cardinality * old_size +
        num_threads * value_type_cardinality * old_size);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        old_cc + (num_threads * value_type_cardinality * reference_cont.size()));
    }
    check_raii_counts();

    // lhs non-empty, rhs non-empty, unequal allocators, no propagation
    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_size = x.size();
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      std::atomic<unsigned> num_transfers{0};

      thread_runner(values, [&x, &values, &num_transfers, &reference_cont](
                              boost::span<value_type> s) {
        (void)s;

        X y(values.size(), hasher(2), key_equal(1), allocator_type(13));
        for (auto const& v : values) {
          y.insert(v);
        }

        BOOST_TEST(
          !boost::allocator_is_always_equal<allocator_type>::type::value);

        BOOST_TEST(!boost::allocator_propagate_on_container_move_assignment<
                   allocator_type>::type::value);

        BOOST_TEST(!y.empty());
        BOOST_TEST(x.get_allocator() != y.get_allocator());

        y = std::move(x);
        if (y.hash_function() == hasher(1)) {
          ++num_transfers;
          test_matches_reference(y, reference_cont);

          BOOST_TEST_EQ(y.key_eq(), key_equal(2));
        } else {
          BOOST_TEST_EQ(y.hash_function(), hasher(2));
          BOOST_TEST_EQ(y.key_eq(), key_equal(1));
        }

        BOOST_TEST(x.empty());

        BOOST_TEST_EQ(x.hash_function(), hasher(2));
        BOOST_TEST_EQ(x.key_eq(), key_equal(1));
        BOOST_TEST(x.get_allocator() != y.get_allocator());
      });

      BOOST_TEST_EQ(num_transfers, 1u);

      BOOST_TEST_EQ(
        raii::destructor,
        2 * value_type_cardinality * old_size + 
        num_threads * value_type_cardinality * old_size);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(
        raii::move_constructor, old_mc + value_type_cardinality * old_size);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        old_cc + (num_threads * value_type_cardinality * reference_cont.size()));
    }
    check_raii_counts();

    // lhs non-empty, rhs non-empty, pocma
    {
      raii::reset_counts();

      pocma_container_type x(
        values.size(), hasher(1), key_equal(2), pocma_allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_size = x.size();
      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      std::atomic<unsigned> num_transfers{0};

      thread_runner(values, [&x, &values, &num_transfers, &reference_cont](
                              boost::span<value_type> s) {
        (void)s;

        pocma_container_type y(
          values.size(), hasher(2), key_equal(1), pocma_allocator_type(13));
        for (auto const& v : values) {
          y.insert(v);
        }

        BOOST_TEST(!y.empty());
        BOOST_TEST(x.get_allocator() != y.get_allocator());

        y = std::move(x);
        if (y.hash_function() == hasher(1)) {
          ++num_transfers;
          test_matches_reference(y, reference_cont);

          BOOST_TEST_EQ(y.key_eq(), key_equal(2));
        } else {
          BOOST_TEST_EQ(y.hash_function(), hasher(2));
          BOOST_TEST_EQ(y.key_eq(), key_equal(1));
        }

        BOOST_TEST(x.empty());

        BOOST_TEST_EQ(x.hash_function(), hasher(2));
        BOOST_TEST_EQ(x.key_eq(), key_equal(1));
        BOOST_TEST(x.get_allocator() == y.get_allocator());
      });

      BOOST_TEST_EQ(num_transfers, 1u);

      BOOST_TEST_EQ(
        raii::destructor,
        value_type_cardinality * old_size +
        num_threads * value_type_cardinality * old_size);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
      BOOST_TEST_EQ(
        raii::copy_constructor,
        old_cc + (num_threads * value_type_cardinality * reference_cont.size()));
    }
    check_raii_counts();

    // self-assign
    {
      raii::reset_counts();

      X x(values.size(), hasher(1), key_equal(2), allocator_type(3));
      for (auto const& v : values) {
        x.insert(v);
      }

      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      thread_runner(
        values, [&x, &reference_cont](boost::span<value_type> s) {
          (void)s;

          x = std::move(x);

          BOOST_TEST(!x.empty());

          BOOST_TEST_EQ(x.hash_function(), hasher(1));
          BOOST_TEST_EQ(x.key_eq(), key_equal(2));
          BOOST_TEST(x.get_allocator() == allocator_type(3));

          test_matches_reference(x, reference_cont);
        });

      BOOST_TEST_EQ(raii::destructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
    }
    check_raii_counts();
  }

  template <class X, class IL>
  void initializer_list_assign(std::pair<X*, IL> p)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto init_list = p.second;
    auto reference_cont = reference_container<X>(
      init_list.begin(), init_list.end());
    auto v = std::vector<value_type>(init_list.begin(), init_list.end());

    {
      raii::reset_counts();
      X x(0, hasher(1), key_equal(2), allocator_type(3));

      thread_runner(v, [&x, &init_list](boost::span<value_type> s) {
        (void)s;
        x = init_list;
      });

      test_matches_reference(x, reference_cont);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(
        raii::copy_constructor, 
        num_threads * value_type_cardinality * x.size());
      BOOST_TEST_EQ(
        raii::destructor,
        (num_threads - 1) * value_type_cardinality * x.size());
      BOOST_TEST_EQ(raii::move_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }
    check_raii_counts();
  }

  template <class X, class GF>
  void initializer_list_assign_gh276(
    X*, GF gen_factory, test::random_generator rg)
  {
    // https://github.com/boostorg/unordered/issues/276

    using replaced_allocator_container = test::replace_allocator<
      X, test::non_default_ctble_allocator<int> >;
    using replaced_allocator_type = 
      typename replaced_allocator_container::allocator_type;
      
    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(4, [&] { return gen(rg); });

    replaced_allocator_container
      x(replaced_allocator_type(0)),
      y(values.begin(), values.end(), replaced_allocator_type(0));

    x = {values[0], values[1], values[2], values[3]};
    BOOST_TEST(x == y);
  }

  template <class X, class GF>
  void insert_and_assign(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();

    std::thread t1, t2, t3;

    boost::compat::latch start_latch(2), end_latch(2);

    auto v1 = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto v2 = v1;
    shuffle_values(v2);

    auto reference_cont = reference_container<X>(v1.begin(), v1.end());

    raii::reset_counts();
    {
      X c1(v1.size(), hasher(1), key_equal(2), allocator_type(3));
      X c2(v2.size(), hasher(1), key_equal(2), allocator_type(3));

      t1 = std::thread([&v1, &c1, &start_latch, &end_latch] {
        start_latch.arrive_and_wait();
        for (auto const& v : v1) {
          c1.insert(v);
        }
        end_latch.arrive_and_wait();
      });

      t2 = std::thread([&v2, &c2, &end_latch, &start_latch] {
        start_latch.arrive_and_wait();
        for (auto const& v : v2) {
          c2.insert(v);
        }
        end_latch.arrive_and_wait();
      });

      std::atomic<unsigned> num_assignments{0};
      t3 = std::thread([&c1, &c2, &end_latch, &num_assignments] {
        while (c1.empty() && c2.empty()) {
          std::this_thread::sleep_for(std::chrono::microseconds(10));
        }

        do {
          c1 = c2;
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          c2 = c1;
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          ++num_assignments;
        } while (!end_latch.try_wait());
      });

      t1.join();
      t2.join();
      t3.join();

      BOOST_TEST_GT(num_assignments, 0u);

      test_fuzzy_matches_reference(c1, reference_cont, rg);
      test_fuzzy_matches_reference(c2, reference_cont, rg);
    }
    check_raii_counts();
  }

  template <class X, class GF>
  void nonconcurrent_move_assign(X*, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    /*
     * basically test that a temporary container is materialized and we
     * move-assign from that
     *
     * we don't need to be super rigorous here because we already have tests for
     * container assignment, we're just testing that a temporary is materialized
     */

    {
      raii::reset_counts();

      nonconcurrent_container<X> nonc(
        values.begin(), values.end(), values.size(),
        hasher(1), key_equal(2), allocator_type(3));

      X x(0, hasher(2), key_equal(1), allocator_type(3));

      BOOST_TEST(nonc.get_allocator() == x.get_allocator());

      x = std::move(nonc);

      BOOST_TEST(nonc.empty());
      BOOST_TEST_EQ(x.size(), reference_cont.size());

      test_fuzzy_matches_reference(x, reference_cont, rg);

      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));

      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(raii::destructor, 0u);
      BOOST_TEST_EQ(raii::move_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }

    check_raii_counts();

    {
      raii::reset_counts();

      X x(values.begin(), values.end(), values.size(), hasher(1),
        key_equal(2), allocator_type(3));

      nonconcurrent_container<X> nonc(
        0, hasher(2), key_equal(1), allocator_type(3));

      BOOST_TEST(nonc.get_allocator() == x.get_allocator());

      nonc = std::move(x);

      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(nonc.size(), reference_cont.size());

      BOOST_TEST_EQ(nonc.hash_function(), hasher(1));
      BOOST_TEST_EQ(nonc.key_eq(), key_equal(2));

      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(raii::destructor, 0u);
      BOOST_TEST_EQ(raii::move_constructor, 0u);
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }

    check_raii_counts();

    {
      raii::reset_counts();

      nonconcurrent_container<X> nonc(
        values.begin(), values.end(), values.size(),
        hasher(1), key_equal(2), allocator_type(3));

      X x(0, hasher(2), key_equal(1), allocator_type(4));

      BOOST_TEST(nonc.get_allocator() != x.get_allocator());

      x = std::move(nonc);

      BOOST_TEST(nonc.empty());
      BOOST_TEST_EQ(x.size(), reference_cont.size());

      test_fuzzy_matches_reference(x, reference_cont, rg);

      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));

      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(
        raii::destructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(
        raii::move_constructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }

    check_raii_counts();

    {
      raii::reset_counts();

      X x(values.begin(), values.end(), values.size(), hasher(1),
        key_equal(2), allocator_type(3));

      nonconcurrent_container<X> nonc(
        0, hasher(2), key_equal(1), allocator_type(4));

      BOOST_TEST(nonc.get_allocator() != x.get_allocator());

      nonc = std::move(x);

      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(nonc.size(), reference_cont.size());

      BOOST_TEST_EQ(nonc.hash_function(), hasher(1));
      BOOST_TEST_EQ(nonc.key_eq(), key_equal(2));

      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(
        raii::destructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(
        raii::move_constructor, value_type_cardinality * reference_cont.size());
      BOOST_TEST_EQ(raii::copy_assignment, 0u);
      BOOST_TEST_EQ(raii::move_assignment, 0u);
    }

    check_raii_counts();
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  copy_assign,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  move_assign,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  initializer_list_assign,
  ((test_map_and_init_list)(test_node_map_and_init_list)
   (test_set_and_init_list)(test_node_set_and_init_list)))

UNORDERED_TEST(
  initializer_list_assign_gh276,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)))

UNORDERED_TEST(
  insert_and_assign,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((init_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  nonconcurrent_move_assign,
  ((test_map)(test_node_map)(test_set)(test_node_set)
   (fancy_test_map)(fancy_test_node_map)(fancy_test_set)(fancy_test_node_set))
  ((init_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
