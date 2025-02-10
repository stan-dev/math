// Copyright (C) 2023 Christian Mazakas
// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "helpers.hpp"

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

test::seed_t initialize_seed{996130204};

using test::default_generator;
using test::limited_range;
using test::sequential;

template <class T> struct pocs_allocator
{
  using propagate_on_container_swap = std::true_type;

  int x_ = -1;

  using value_type = T;

  pocs_allocator() = default;
  pocs_allocator(pocs_allocator const&) = default;
  pocs_allocator(pocs_allocator&&) = default;

  pocs_allocator(int const x) : x_{x} {}

  pocs_allocator& operator=(pocs_allocator const& rhs)
  {
    if (this != &rhs) {
      x_ = rhs.x_;
    }
    return *this;
  }

  template <class U> pocs_allocator(pocs_allocator<U> const& rhs) : x_{rhs.x_}
  {
  }

  T* allocate(std::size_t n)
  {
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) { ::operator delete(p); }

  bool operator==(pocs_allocator const& rhs) const { return x_ == rhs.x_; }
  bool operator!=(pocs_allocator const& rhs) const { return x_ != rhs.x_; }

  friend void swap(pocs_allocator& lhs, pocs_allocator& rhs) noexcept
  {
    std::swap(lhs.x_, rhs.x_);
  }
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

template <class T> struct is_nothrow_member_swappable
{
  static bool const value =
    noexcept(std::declval<T&>().swap(std::declval<T&>()));
};

BOOST_STATIC_ASSERT(is_nothrow_member_swappable<
  replace_allocator<map_type, std::allocator> >::value);

BOOST_STATIC_ASSERT(is_nothrow_member_swappable<
  replace_allocator<map_type, pocs_allocator> >::value);

BOOST_STATIC_ASSERT(!is_nothrow_member_swappable<map_type>::value);

BOOST_STATIC_ASSERT(is_nothrow_member_swappable<
  replace_allocator<set_type, std::allocator> >::value);

BOOST_STATIC_ASSERT(is_nothrow_member_swappable<
  replace_allocator<set_type, pocs_allocator> >::value);

BOOST_STATIC_ASSERT(!is_nothrow_member_swappable<set_type>::value);

namespace {
  struct
  {
    template <class T> void operator()(T& x1, T& x2) const { x1.swap(x2); }
  } member_fn_swap;

  struct
  {
    template <class T> void operator()(T& x1, T& x2) const
    {
      using boost::unordered::swap;
      swap(x1, x2);
    }
  } free_fn_swap;

  template <class X, class F, class GF>
  void swap_tests(X*, F swapper, GF gen_factory, test::random_generator rg)
  {
    using value_type = typename X::value_type;
    using allocator_type = typename X::allocator_type;

    bool const pocs =
      boost::allocator_propagate_on_container_swap<
        allocator_type>::type::value;

    auto gen = gen_factory.template get<X>();
    auto vals1 = make_random_values(1024 * 8, [&] { return gen(rg); });
    auto vals2 = make_random_values(1024 * 4, [&] { return gen(rg); });

    auto reference_cont1 = reference_container<X>(vals1.begin(), vals1.end());
    auto reference_cont2 = reference_container<X>(vals2.begin(), vals2.end());

    {
      raii::reset_counts();

      X x1(vals1.begin(), vals1.end(), vals1.size(), hasher(1), key_equal(2),
        allocator_type(3));

      X x2(vals2.begin(), vals2.end(), vals2.size(), hasher(2), key_equal(1),
        pocs ? allocator_type(4) : allocator_type(3));

      if (pocs) {
        BOOST_TEST(x1.get_allocator() != x2.get_allocator());
      } else {
        BOOST_TEST(x1.get_allocator() == x2.get_allocator());
      }

      auto const old_cc = +raii::copy_constructor;
      auto const old_mc = +raii::move_constructor;

      thread_runner(vals1, [&x1, &x2, swapper](boost::span<value_type> s) {
        (void)s;

        swapper(x1, x2);
        swapper(x2, x1);
      });

      BOOST_TEST_EQ(raii::copy_constructor, old_cc);
      BOOST_TEST_EQ(raii::move_constructor, old_mc);

      if (pocs) {
        if (x1.get_allocator() == allocator_type(3)) {
          BOOST_TEST(x2.get_allocator() == allocator_type(4));
        } else {
          BOOST_TEST(x1.get_allocator() == allocator_type(4));
          BOOST_TEST(x2.get_allocator() == allocator_type(3));
        }
      } else {
        BOOST_TEST(x1.get_allocator() == allocator_type(3));
        BOOST_TEST(x1.get_allocator() == x2.get_allocator());
      }

      if (x1.size() == reference_cont1.size()) {
        test_matches_reference(x1, reference_cont1);
        test_matches_reference(x2, reference_cont2);

        BOOST_TEST_EQ(x1.hash_function(), hasher(1));
        BOOST_TEST_EQ(x1.key_eq(), key_equal(2));

        BOOST_TEST_EQ(x2.hash_function(), hasher(2));
        BOOST_TEST_EQ(x2.key_eq(), key_equal(1));
      } else {
        test_matches_reference(x2, reference_cont1);
        test_matches_reference(x1, reference_cont2);

        BOOST_TEST_EQ(x1.hash_function(), hasher(2));
        BOOST_TEST_EQ(x1.key_eq(), key_equal(1));

        BOOST_TEST_EQ(x2.hash_function(), hasher(1));
        BOOST_TEST_EQ(x2.key_eq(), key_equal(2));
      }
    }
    check_raii_counts();
  }

  template <class X, class F, class GF>
  void insert_and_swap(
    X*, F swapper, GF gen_factory, test::random_generator rg)
  {
    using allocator_type = typename X::allocator_type;

    auto gen = gen_factory.template get<X>();
    auto vals1 = make_random_values(1024 * 8, [&] { return gen(rg); });
    auto vals2 = make_random_values(1024 * 4, [&] { return gen(rg); });

    {
      raii::reset_counts();

      X x1(vals1.size(), hasher(1), key_equal(2), allocator_type(3));
      X x2(vals2.size(), hasher(2), key_equal(1), allocator_type(3));

      std::thread t1, t2, t3;
      boost::compat::latch l(2);

      std::mutex m;
      std::condition_variable cv;
      std::atomic_bool done1{false}, done2{false};
      std::atomic<unsigned> num_swaps{0};
      bool ready = false;

      t1 = std::thread([&x1, &vals1, &l, &done1, &cv, &ready, &m] {
        l.arrive_and_wait();

        for (std::size_t idx = 0; idx < vals1.size(); ++idx) {
          auto const& val = vals1[idx];
          x1.insert(val);
          if (idx % (vals1.size() / 128) == 0) {
            {
              std::unique_lock<std::mutex> lk(m);
              ready = true;
            }
            cv.notify_all();
          }
          std::this_thread::yield();
        }

        done1 = true;
        {
          std::unique_lock<std::mutex> lk(m);
          ready = true;
        }
        cv.notify_all();
      });

      t2 = std::thread([&x2, &vals2, &l, &done2, &ready, &cv, &m] {
        l.arrive_and_wait();

        for (std::size_t idx = 0; idx < vals2.size(); ++idx) {
          auto const& val = vals2[idx];
          x2.insert(val);
          if (idx % 100 == 0) {
            std::this_thread::yield();
          }
        }

        done2 = true;
        {
          std::unique_lock<std::mutex> lk(m);
          ready = true;
        }
        cv.notify_all();
      });

      t3 = std::thread(
        [&x1, &x2, &m, &cv, &done1, &done2, &num_swaps, swapper, &ready] {
          do {
            {
              std::unique_lock<std::mutex> lk(m);
              cv.wait(lk, [&ready] { return ready; });
              ready = false;
            }
            swapper(x1, x2);
            ++num_swaps;
            std::this_thread::yield();
          } while (!done1 || !done2);

          BOOST_TEST(done1);
          BOOST_TEST(done2);
        });

      t1.join();
      t2.join();
      t3.join();

      BOOST_TEST_GT(num_swaps, 0u);

      if (x1.hash_function() == hasher(1)) {
        BOOST_TEST_EQ(x1.key_eq(), key_equal(2));

        BOOST_TEST_EQ(x2.hash_function(), hasher(2));
        BOOST_TEST_EQ(x2.key_eq(), key_equal(1));
      } else {
        BOOST_TEST_EQ(x1.hash_function(), hasher(2));
        BOOST_TEST_EQ(x1.key_eq(), key_equal(1));

        BOOST_TEST_EQ(x2.hash_function(), hasher(1));
        BOOST_TEST_EQ(x2.key_eq(), key_equal(2));
      }
    }

    check_raii_counts();
  }

  map_type* map;
  replace_allocator<map_type, pocs_allocator>* pocs_map;

  node_map_type* node_map;
  replace_allocator<node_map_type, pocs_allocator>* pocs_node_map;

  set_type* set;
  replace_allocator<set_type, pocs_allocator>* pocs_set;

  node_set_type* node_set;
  replace_allocator<node_set_type, pocs_allocator>* pocs_node_set;

} // namespace

// clang-format off
UNORDERED_TEST(
  swap_tests,
  ((map)(pocs_map)(node_map)(pocs_node_map)
   (set)(pocs_set)(node_set)(pocs_node_set))
  ((member_fn_swap)(free_fn_swap))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))

UNORDERED_TEST(insert_and_swap,
  ((map)(node_map)(set)(node_set))
  ((member_fn_swap)(free_fn_swap))
  ((value_type_generator_factory))
  ((default_generator)(sequential)(limited_range)))
// clang-format on

RUN_TESTS()
