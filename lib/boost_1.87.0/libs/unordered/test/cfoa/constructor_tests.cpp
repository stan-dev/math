// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

test::seed_t initialize_seed(4122023);

using test::default_generator;
using test::limited_range;
using test::sequential;

template <class T> struct soccc_allocator
{
  int x_ = -1;

  using value_type = T;

  soccc_allocator() = default;
  soccc_allocator(soccc_allocator const&) = default;
  soccc_allocator(soccc_allocator&&) = default;

  soccc_allocator(int const x) : x_{x} {}

  template <class U> soccc_allocator(soccc_allocator<U> const& rhs) : x_{rhs.x_}
  {
  }

  T* allocate(std::size_t n)
  {
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) { ::operator delete(p); }

  soccc_allocator select_on_container_copy_construction() const
  {
    return {x_ + 1};
  }

  bool operator==(soccc_allocator const& rhs) const { return x_ == rhs.x_; }
  bool operator!=(soccc_allocator const& rhs) const { return x_ != rhs.x_; }
};

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

map_type* test_map;
node_map_type* test_node_map;
set_type* test_set;
node_set_type* test_node_set;

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

namespace {
  template <class X>
  void default_constructor(X*)
  {
    X x;
    BOOST_TEST(x.empty());
    BOOST_TEST_EQ(x.size(), 0u);
  }

  template <class X>
  void bucket_count_with_hasher_key_equal_and_allocator(X*)
  {
    using allocator_type = typename X::allocator_type;

    raii::reset_counts();
    {
      X x(0);

      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
    }

    {
      X x(0, hasher(1));

      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal());
    }

    {
      X x(0, hasher(1), key_equal(2));

      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
    }

    {
      X x(0, hasher(1), key_equal(2), allocator_type{});

      BOOST_TEST(x.empty());
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type{});
    }
  }

  template <class X>
  void soccc(X*)
  {
    raii::reset_counts();

    replace_allocator<X, soccc_allocator> x, y(x);

    BOOST_TEST_EQ(y.hash_function(), x.hash_function());
    BOOST_TEST_EQ(y.key_eq(), x.key_eq());
    BOOST_TEST(y.get_allocator() != x.get_allocator());
  }

  template <class X, class GF>
  void from_iterator_range(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());
    raii::reset_counts();

    {
      X x(values.begin(), values.end());

      test_matches_reference(x, reference_cont);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      X x(values.begin(), values.end(), 0);

      test_matches_reference(x, reference_cont);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      X x(values.begin(), values.end(), 0, hasher(1));

      test_matches_reference(x, reference_cont);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      X x(values.begin(), values.end(), 0, hasher(1), key_equal(2));

      test_matches_reference(x, reference_cont);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    {
      X x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{});

      test_matches_reference(x, reference_cont);
      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type{});
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }
    }

    check_raii_counts();
  }

  template <class X, class GF>
  void copy_constructor(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    {
      X x(0, hasher(1), key_equal(2), allocator_type{});
      X y(x);

      BOOST_TEST_EQ(y.size(), x.size());
      BOOST_TEST_EQ(y.hash_function(), x.hash_function());
      BOOST_TEST_EQ(y.key_eq(), x.key_eq());
      BOOST_TEST(y.get_allocator() == x.get_allocator());
    }

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      X x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{});

      thread_runner(
        values, [&x, &reference_cont](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;
          X y(x);

          test_matches_reference(x, reference_cont);
          test_matches_reference(y, reference_cont);
          BOOST_TEST_EQ(y.size(), x.size());
          BOOST_TEST_EQ(y.hash_function(), x.hash_function());
          BOOST_TEST_EQ(y.key_eq(), x.key_eq());
          BOOST_TEST(y.get_allocator() == x.get_allocator());
        });
    }

    check_raii_counts();

    raii::reset_counts();

    {
      allocator_type a;

      X x(values.begin(), values.end(), 0, hasher(1), key_equal(2), a);

      thread_runner(
        values, [&x, &reference_cont, a](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;
          X y(x, a);

          test_matches_reference(x, reference_cont);
          test_matches_reference(y, reference_cont);
          BOOST_TEST_EQ(y.size(), x.size());
          BOOST_TEST_EQ(y.hash_function(), x.hash_function());
          BOOST_TEST_EQ(y.key_eq(), x.key_eq());
          BOOST_TEST(y.get_allocator() == x.get_allocator());
        });
    }

    check_raii_counts();
  }

  template <class X, class GF>
  void copy_constructor_with_insertion(X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());
    raii::reset_counts();

    std::mutex m;
    std::condition_variable cv;
    bool ready = false;

    {
      X x(0, hasher(1), key_equal(2), allocator_type{});

      auto f = [&x, &values, &m, &cv, &ready] {
        {
          std::lock_guard<std::mutex> guard(m);
          ready = true;
        }
        cv.notify_all();

        for (auto const& val : values) {
          x.insert(val);
        }
      };

      std::thread t1(f);
      std::thread t2(f);

      thread_runner(
        values, [&x, &reference_cont, &values, rg, &m, &cv, &ready](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [&] { return ready; });
          }

          X y(x);

          BOOST_TEST_LE(y.size(), values.size());
          BOOST_TEST_EQ(y.hash_function(), x.hash_function());
          BOOST_TEST_EQ(y.key_eq(), x.key_eq());
          BOOST_TEST(y.get_allocator() == x.get_allocator());

          x.visit_all([&reference_cont, rg](
                        typename X::value_type const& val) {
            BOOST_TEST(reference_cont.contains(get_key(val)));
            if (rg == sequential) {
              BOOST_TEST_EQ(val, *reference_cont.find(get_key(val)));
            }
          });
        });

      t1.join();
      t2.join();
    }

    check_raii_counts();
  }

  template <class X, class GF>
  void move_constructor(X*, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    using allocator_type = typename X::allocator_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;

    {
      X x(0, hasher(1), key_equal(2), allocator_type{});
      auto const old_size = x.size();

      X y(std::move(x));

      BOOST_TEST_EQ(y.size(), old_size);
      BOOST_TEST_EQ(y.hash_function(), hasher(1));
      BOOST_TEST_EQ(y.key_eq(), key_equal(2));

      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());

      BOOST_TEST(y.get_allocator() == x.get_allocator());
    }

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      X x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{});

      std::atomic_uint num_transfers{0};

      auto const old_mc = +raii::move_constructor;

      thread_runner(
        values, [&x, &reference_cont, &num_transfers](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          auto const old_size = x.size();
          X y(std::move(x));

          if (!y.empty()) {
            ++num_transfers;

            test_matches_reference(y, reference_cont);
            BOOST_TEST_EQ(y.size(), old_size);
            BOOST_TEST_EQ(y.hash_function(), hasher(1));
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.size(), 0u);
            BOOST_TEST_EQ(y.hash_function(), hasher());
            BOOST_TEST_EQ(y.key_eq(), key_equal());
          }

          BOOST_TEST_EQ(x.size(), 0u);
          BOOST_TEST_EQ(x.hash_function(), hasher());
          BOOST_TEST_EQ(x.key_eq(), key_equal());

          BOOST_TEST(y.get_allocator() == x.get_allocator());
        });

      BOOST_TEST_EQ(num_transfers, 1u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
    }

    check_raii_counts();

    // allocator-aware move constructor, unequal allocators
    raii::reset_counts();

    {
      X x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{1});

      std::atomic_uint num_transfers{0};

      auto const old_mc = +raii::move_constructor;
      auto const old_size = x.size();

      thread_runner(
        values, [&x, &reference_cont, &num_transfers, old_size](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          auto a = allocator_type{2};
          BOOST_TEST(a != x.get_allocator());

          X y(std::move(x), a);

          if (!y.empty()) {
            ++num_transfers;

            test_matches_reference(y, reference_cont);
            BOOST_TEST_EQ(y.size(), old_size);
            BOOST_TEST_EQ(y.hash_function(), hasher(1));
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.size(), 0u);
            BOOST_TEST_EQ(y.hash_function(), hasher());
            BOOST_TEST_EQ(y.key_eq(), key_equal());
          }

          BOOST_TEST_EQ(x.size(), 0u);
          BOOST_TEST_EQ(x.hash_function(), hasher());
          BOOST_TEST_EQ(x.key_eq(), key_equal());

          BOOST_TEST(y.get_allocator() != x.get_allocator());
          BOOST_TEST(y.get_allocator() == a);
        });

      BOOST_TEST_EQ(num_transfers, 1u);
      BOOST_TEST_EQ(
        raii::move_constructor, old_mc + (value_type_cardinality * old_size));
    }

    check_raii_counts();

    // allocator-aware move constructor, equal allocators
    raii::reset_counts();

    {
      X x(values.begin(), values.end(), 0, hasher(1), key_equal(2),
        allocator_type{1});

      std::atomic_uint num_transfers{0};

      auto const old_mc = +raii::move_constructor;
      auto const old_size = x.size();

      thread_runner(
        values, [&x, &reference_cont, &num_transfers, old_size](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          auto a = allocator_type{1};
          BOOST_TEST(a == x.get_allocator());

          X y(std::move(x), a);

          if (!y.empty()) {
            ++num_transfers;

            test_matches_reference(y, reference_cont);
            BOOST_TEST_EQ(y.size(), old_size);
            BOOST_TEST_EQ(y.hash_function(), hasher(1));
            BOOST_TEST_EQ(y.key_eq(), key_equal(2));
          } else {
            BOOST_TEST_EQ(y.size(), 0u);
            BOOST_TEST_EQ(y.hash_function(), hasher());
            BOOST_TEST_EQ(y.key_eq(), key_equal());
          }

          BOOST_TEST_EQ(x.size(), 0u);
          BOOST_TEST_EQ(x.hash_function(), hasher());
          BOOST_TEST_EQ(x.key_eq(), key_equal());

          BOOST_TEST(y.get_allocator() == x.get_allocator());
          BOOST_TEST(y.get_allocator() == a);
        });

      BOOST_TEST_EQ(num_transfers, 1u);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);
    }

    check_raii_counts();
  }

  template <class X, class GF>
  void move_constructor_with_insertion(
    X*, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    std::mutex m;
    std::condition_variable cv;
    bool ready = false;

    {
      X x(0, hasher(1), key_equal(2), allocator_type{});

      std::atomic_uint num_transfers{0};

      std::thread t1([&x, &values] {
        for (auto const& val : values) {
          x.insert(val);
        }
      });

      std::thread t2([&x, &m, &cv, &ready] {
        while (x.empty()) {
          std::this_thread::yield();
        }

        {
          std::lock_guard<std::mutex> guard(m);
          ready = true;
        }
        cv.notify_all();
      });

      thread_runner(
        values, [&x, &reference_cont, &num_transfers, rg, &m, &ready, &cv](
                  boost::span<span_value_type<decltype(values)> > s) {
          (void)s;

          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [&] { return ready; });
          }

          X y(std::move(x));

          if (!y.empty()) {
            ++num_transfers;
            y.cvisit_all([&reference_cont, rg](value_type const& val) {
              BOOST_TEST(reference_cont.contains(get_key(val)));
              if (rg == sequential) {
                BOOST_TEST_EQ(
                  val, *reference_cont.find(get_key(val)));
              }
            });
          }
        });

      t1.join();
      t2.join();

      BOOST_TEST_GE(num_transfers, 1u);
    }

    check_raii_counts();
  }

  template <class X, class GF>
  void iterator_range_with_allocator(
    X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      allocator_type a;
      X x(values.begin(), values.end(), a);

      BOOST_TEST_GT(x.size(), 0u);
      BOOST_TEST_LE(x.size(), values.size());
      if (rg == sequential) {
        BOOST_TEST_EQ(x.size(), values.size());
      }

      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());

      BOOST_TEST(x.get_allocator() == a);

      test_fuzzy_matches_reference(x, reference_cont, rg);
    }

    check_raii_counts();
  }

  template <class X>
  void explicit_allocator(X*)
  {
    using allocator_type = typename X::allocator_type;

    raii::reset_counts();

    {
      allocator_type a;
      X x(a);

      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());

      BOOST_TEST(x.get_allocator() == a);
    }
  }

  template <class X, class IL>
  void initializer_list_with_all_params(std::pair<X*, IL> p)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality = 
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto init_list = p.second;

    {
      raii::reset_counts();

      X x(init_list, 0, hasher(1), key_equal(2), allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * init_list.size() / 2u);
      BOOST_TEST_EQ(
        raii::move_constructor, 0u);
    }
    check_raii_counts();

    {
      raii::reset_counts();

      X x(init_list, allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * init_list.size() / 2u);
      BOOST_TEST_EQ(
        raii::move_constructor, 0u);
    }
    check_raii_counts();

    {
      raii::reset_counts();

      X x(init_list, 0, allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * init_list.size() / 2u);
      BOOST_TEST_EQ(
        raii::move_constructor, 0u);
    }
    check_raii_counts();

    {
      raii::reset_counts();

      X x(init_list, 0, hasher(1), allocator_type(3));

      BOOST_TEST_EQ(x.size(), 11u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST_EQ(raii::default_constructor, 0u);
      BOOST_TEST_EQ(
        raii::copy_constructor, value_type_cardinality * init_list.size() / 2u);
      BOOST_TEST_EQ(
        raii::move_constructor, 0u);
    }
    check_raii_counts();
  }

  template <class X>
  void bucket_count_and_allocator(X*)
  {
    using allocator_type = typename X::allocator_type;

    raii::reset_counts();

    {
      X x(0, allocator_type(3));
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));
    }

    {
      X x(4096, allocator_type(3));
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));
    }
  }

  template <class X>
  void bucket_count_with_hasher_and_allocator(X*)
  {
    using allocator_type = typename X::allocator_type;

    raii::reset_counts();

    {
      X x(0, hasher(1), allocator_type(3));
      BOOST_TEST_EQ(x.size(), 0u);
      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == allocator_type(3));
    }
  }

  template <class X, class GF>
  void iterator_range_with_bucket_count_and_allocator(
    X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      allocator_type a(3);
      X x(values.begin(), values.end(), 0, a);
      test_fuzzy_matches_reference(x, reference_cont, rg);

      BOOST_TEST_EQ(x.hash_function(), hasher());
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == a);
    }

    check_raii_counts();
  }

  template <class X, class GF>
  void iterator_range_with_bucket_count_hasher_and_allocator(
    X*, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      allocator_type a(3);
      hasher hf(1);
      X x(values.begin(), values.end(), 0, hf, a);
      test_fuzzy_matches_reference(x, reference_cont, rg);

      BOOST_TEST_EQ(x.hash_function(), hf);
      BOOST_TEST_EQ(x.key_eq(), key_equal());
      BOOST_TEST(x.get_allocator() == a);
    }

    check_raii_counts();
  }

  template <class X, class GF>
  void nonconcurrent_constructor(X*, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    static constexpr auto value_type_cardinality =
      value_cardinality<value_type>::value;
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto values = make_random_values(1024 * 16, [&] { return gen(rg); });
    auto reference_cont = reference_container<X>(values.begin(), values.end());
    auto reference_nonc = 
      nonconcurrent_container<X>(values.begin(), values.end());

    raii::reset_counts();

    {
      nonconcurrent_container<X> nonc(
        values.begin(), values.end(), reference_cont.size(), hasher(1),
        key_equal(2), allocator_type(3));

      auto const old_dc = +raii::default_constructor;
      auto const old_mc = +raii::move_constructor;
      auto const old_cc = +raii::copy_constructor;

      BOOST_TEST_EQ(old_dc, 0u);
      BOOST_TEST_EQ(old_mc, 0u);
      BOOST_TEST_EQ(old_cc, value_type_cardinality * nonc.size());

      X x(std::move(nonc));

      test_fuzzy_matches_reference(x, reference_cont, rg);

      BOOST_TEST_EQ(+raii::default_constructor, old_dc);
      BOOST_TEST_EQ(+raii::move_constructor, old_mc);
      BOOST_TEST_EQ(+raii::copy_constructor, old_cc);

      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST(nonc.empty());
    }

    check_raii_counts();

    {
      nonconcurrent_container<X> nonc(
        0, hasher(1), key_equal(2), allocator_type(3));

      X x(std::move(nonc));

      BOOST_TEST(x.empty());

      BOOST_TEST_EQ(x.hash_function(), hasher(1));
      BOOST_TEST_EQ(x.key_eq(), key_equal(2));
      BOOST_TEST(x.get_allocator() == allocator_type(3));

      BOOST_TEST(nonc.empty());
    }

    check_raii_counts();

    {
      X x(values.begin(), values.end(), reference_cont.size(),
        hasher(1), key_equal(2), allocator_type(3));

      auto const old_dc = +raii::default_constructor;
      auto const old_mc = +raii::move_constructor;
      auto const old_cc = +raii::copy_constructor;

      BOOST_TEST_EQ(old_dc, 0u);
      BOOST_TEST_EQ(old_mc, 0u);
      BOOST_TEST_EQ(old_cc, 2u * value_type_cardinality * x.size());

      nonconcurrent_container<X> nonc(std::move(x));

      BOOST_TEST(nonc == reference_nonc);

      BOOST_TEST_EQ(+raii::default_constructor, old_dc);
      BOOST_TEST_EQ(+raii::move_constructor, old_mc);
      BOOST_TEST_EQ(+raii::copy_constructor, old_cc);

      BOOST_TEST_EQ(nonc.hash_function(), hasher(1));
      BOOST_TEST_EQ(nonc.key_eq(), key_equal(2));
      BOOST_TEST(nonc.get_allocator() == allocator_type(3));

      BOOST_TEST(x.empty());
    }

    check_raii_counts();

    {
      X x(0, hasher(1), key_equal(2), allocator_type(3));

      nonconcurrent_container<X> nonc(std::move(x));

      BOOST_TEST(nonc.empty());

      BOOST_TEST_EQ(nonc.hash_function(), hasher(1));
      BOOST_TEST_EQ(nonc.key_eq(), key_equal(2));
      BOOST_TEST(nonc.get_allocator() == allocator_type(3));

      BOOST_TEST(x.empty());
    }

    check_raii_counts();
  }

} // namespace

// clang-format off
UNORDERED_TEST(
  default_constructor,
  ((test_map)(test_node_map)(test_set)(test_node_set)))

UNORDERED_TEST(
  bucket_count_with_hasher_key_equal_and_allocator,
  ((test_map)(test_node_map)(test_set)(test_node_set)))

UNORDERED_TEST(
  soccc,
  ((test_map)(test_node_map)(test_set)(test_node_set)))

UNORDERED_TEST(
  from_iterator_range,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  copy_constructor,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  copy_constructor_with_insertion,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  move_constructor,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  move_constructor_with_insertion,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  iterator_range_with_allocator,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  explicit_allocator,
  ((test_map)(test_node_map)(test_set)(test_node_set)))

UNORDERED_TEST(
  initializer_list_with_all_params,
  ((test_map_and_init_list)(test_node_map_and_init_list)
   (test_set_and_init_list)(test_node_set_and_init_list)))

UNORDERED_TEST(
  bucket_count_and_allocator,
  ((test_map)(test_node_map)(test_set)(test_node_set)))

UNORDERED_TEST(
  bucket_count_with_hasher_and_allocator,
  ((test_map)(test_node_map)(test_set)(test_node_set)))
  
UNORDERED_TEST(
  iterator_range_with_bucket_count_and_allocator,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  iterator_range_with_bucket_count_hasher_and_allocator,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(
  nonconcurrent_constructor,
  ((test_map)(test_node_map)(test_set)(test_node_set))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
