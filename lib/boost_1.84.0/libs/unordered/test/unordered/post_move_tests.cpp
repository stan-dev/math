
// Copyright (C) 2022-2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or move at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/equivalent.hpp"
#include "../helpers/invariants.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/test.hpp"
#include "../helpers/tracker.hpp"
#include "../helpers/unordered.hpp"
#include "../objects/cxx11_allocator.hpp"
#include "../objects/test.hpp"

#include <boost/core/ignore_unused.hpp>
#include <iterator>

#if defined(BOOST_MSVC)
#pragma warning(disable : 4127) // conditional expression is constant
#endif

namespace move_tests {
  test::seed_t initialize_seed(98624);
#define BOOST_UNORDERED_TEST_MOVING 1

  template <class T> T empty(T*) { return T(); }

  template <class T>
  T create(test::random_values<T> const& v, test::object_count& count)
  {
    T x(v.begin(), v.end());
    count = test::global_object_count;
    return x;
  }

  template <class T>
  T create(test::random_values<T> const& v, test::object_count& count,
    typename T::hasher hf, typename T::key_equal eq,
    typename T::allocator_type al, float mlf)
  {
    T x(0, hf, eq, al);
    x.max_load_factor(mlf);
    x.insert(v.begin(), v.end());
    count = test::global_object_count;
    return x;
  }

  template <class T> T const& get_key(T const& t) { return t; }

  template <class K, class V> K const& get_key(std::pair<K const, V> const& kv)
  {
    return kv.first;
  }

  template <class T> T const& get_value(T const& t) { return t; }

  template <class K, class V>
  K const& get_value(std::pair<K const, V> const& kv)
  {
    return kv.second;
  }

  template <class T>
  static void insert_range(T& y, test::random_values<T> const& v)
  {
    y.insert(v.begin(), v.end());
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T>
  static void insert_single(T& y, test::random_values<T> const& v)
  {
    y.insert(*v.begin());
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T>
  static void insert_single_hint(T& y, test::random_values<T> const& v)
  {
    y.insert(y.end(), *v.begin());
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T> struct insert_or_assign_invoker
  {
    void operator()(T&, test::random_values<T> const&) {}
  };

  template <class Key, class T, class Hash, class KeyEqual, class Allocator>
  struct insert_or_assign_invoker<
    boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> >
  {
    void operator()(boost::unordered_map<Key, T, Hash, KeyEqual, Allocator>& y,
      test::random_values<
        boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> > const& v)
    {
      typedef typename boost::unordered_map<Key, T, Hash, KeyEqual,
        Allocator>::size_type size_type;

      y.insert_or_assign(get_key(*v.begin()), get_value(*v.begin()));
      BOOST_TEST_EQ(
        y.size(), static_cast<size_type>(std::distance(y.begin(), y.end())));
    }
  };

  template <class T>
  static void insert_or_assign(T& y, test::random_values<T> const& v)
  {
    insert_or_assign_invoker<T>()(y, v);
  }

  template <class T> struct insert_or_assign_hint_invoker
  {
    void operator()(T&, test::random_values<T> const&) {}
  };

  template <class Key, class T, class Hash, class KeyEqual, class Allocator>
  struct insert_or_assign_hint_invoker<
    boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> >
  {
    void operator()(boost::unordered_map<Key, T, Hash, KeyEqual, Allocator>& y,
      test::random_values<
        boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> > const& v)
    {
      typedef typename boost::unordered_map<Key, T, Hash, KeyEqual,
        Allocator>::size_type size_type;
      y.insert_or_assign(y.end(), get_key(*v.begin()), get_value(*v.begin()));
      BOOST_TEST_EQ(
        y.size(), static_cast<size_type>(std::distance(y.begin(), y.end())));
    }
  };

  template <class T>
  static void insert_or_assign_hint(T& y, test::random_values<T> const& v)
  {
    insert_or_assign_hint_invoker<T>()(y, v);
  }

  template <class T> struct try_emplace_invoker
  {
    void operator()(T&, test::random_values<T> const&) {}
  };

  template <class Key, class T, class Hash, class KeyEqual, class Allocator>
  struct try_emplace_invoker<
    boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> >
  {
    void operator()(boost::unordered_map<Key, T, Hash, KeyEqual, Allocator>& y,
      test::random_values<
        boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> > const& v)
    {
      typedef typename boost::unordered_map<Key, T, Hash, KeyEqual,
        Allocator>::size_type size_type;
      y.try_emplace(get_key(*v.begin()), get_value(*v.begin()));
      BOOST_TEST_EQ(
        y.size(), static_cast<size_type>(std::distance(y.begin(), y.end())));
    }
  };

  template <class T>
  static void try_emplace(T& y, test::random_values<T> const& v)
  {
    try_emplace_invoker<T>()(y, v);
  }

  template <class T> struct try_emplace_hint_invoker
  {
    void operator()(T&, test::random_values<T> const&) {}
  };

  template <class Key, class T, class Hash, class KeyEqual, class Allocator>
  struct try_emplace_hint_invoker<
    boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> >
  {
    void operator()(boost::unordered_map<Key, T, Hash, KeyEqual, Allocator>& y,
      test::random_values<
        boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> > const& v)
    {
      typedef typename boost::unordered_map<Key, T, Hash, KeyEqual,
        Allocator>::size_type size_type;
      y.try_emplace(y.end(), get_key(*v.begin()), get_value(*v.begin()));
      BOOST_TEST_EQ(
        y.size(), static_cast<size_type>(std::distance(y.begin(), y.end())));
    }
  };

  template <class T>
  static void try_emplace_hint(T& y, test::random_values<T> const& v)
  {
    try_emplace_hint_invoker<T>()(y, v);
  }

  template <class T> struct at_invoker
  {
    void operator()(T&, test::random_values<T> const&) {}
  };

  template <class Key, class T, class Hash, class KeyEqual, class Allocator>
  struct at_invoker<boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> >
  {
    void operator()(boost::unordered_map<Key, T, Hash, KeyEqual, Allocator>& y,
      test::random_values<
        boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> > const& v)
    {
      BOOST_TRY { y.at(get_key(*v.begin())); }
      BOOST_CATCH(...) {}
      BOOST_CATCH_END
    }
  };

  template <class T> static void at(T& y, test::random_values<T> const& v)
  {
    at_invoker<T>()(y, v);
  }

  template <class T> struct index_operator_invoker
  {
    void operator()(T&, test::random_values<T> const&) {}
  };

  template <class Key, class T, class Hash, class KeyEqual, class Allocator>
  struct index_operator_invoker<
    boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> >
  {
    void operator()(boost::unordered_map<Key, T, Hash, KeyEqual, Allocator>& y,
      test::random_values<
        boost::unordered_map<Key, T, Hash, KeyEqual, Allocator> > const& v)
    {
      typedef typename boost::unordered_map<Key, T, Hash, KeyEqual,
        Allocator>::size_type size_type;
      y[get_key(*v.begin())] = get_value(*v.begin());
      BOOST_TEST_EQ(
        y.size(), static_cast<size_type>(std::distance(y.begin(), y.end())));
    }
  };

  template <class T>
  static void index_operator(T& y, test::random_values<T> const& v)
  {
    index_operator_invoker<T>()(y, v);
  }

  template <class T> static void clear(T& y, test::random_values<T> const&)
  {
    y.clear();
    BOOST_TEST(y.empty());
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T> static void capacity(T& y, test::random_values<T> const&)
  {
    (void)y.empty();
    (void)y.size();
    (void)y.max_size();
    (void)y.load_factor();
    (void)y.max_load_factor();
    (void)y.hash_function();
    (void)y.key_eq();
    (void)y.get_allocator();
  }

  template <class T> static void iterators(T& y, test::random_values<T> const&)
  {
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T>
  static void erase_range(T& y, test::random_values<T> const&)
  {
    y.erase(y.begin(), y.end());
    BOOST_TEST(y.empty());
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T>
  static void erase_key(T& y, test::random_values<T> const& v)
  {
    y.erase(get_key(*v.begin()));
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T> static void lookup(T& y, test::random_values<T> const& v)
  {
    (void)y.count(get_key(*v.begin()));
    (void)y.find(get_key(*v.begin()));
    (void)y.contains(get_key(*v.begin()));
    (void)y.equal_range(get_key(*v.begin()));
  }

  template <class T> static void reserve(T& y, test::random_values<T> const&)
  {
    y.reserve(1337);
    BOOST_TEST_GT(y.bucket_count(), 1337u);
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T>
  static void copy_assignment(T& y, test::random_values<T> const& v)
  {
    T x(v.begin(), v.end());
    y = x;
    BOOST_TEST_EQ(y.size(), x.size());
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T>
  static void move_assignment(T& y, test::random_values<T> const& v)
  {
    T x(v.begin(), v.end());
    std::size_t const size = x.size();
    y = std::move(x);
    BOOST_TEST_GE(y.size(), size);
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T> static void equal(T& y, test::random_values<T> const& v)
  {
    T x(v.begin(), v.end());
    (void)(y == x);
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

#ifndef BOOST_UNORDERED_FOA_TESTS
  template <class T> static void extract(T& y, test::random_values<T> const& v)
  {
    (void)y.extract(get_key(*v.begin()));
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }
#endif

  template <class T> static void merge(T& y, test::random_values<T> const& v)
  {
    T x(v.begin(), v.end());
    if (y.get_allocator() == x.get_allocator()) {
      y.merge(x);
    }
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class X> bool pred(X const&) { return true; }

  template <class T>
  static void erase_with_pred(T& y, test::random_values<T> const&)
  {
    erase_if(y, pred<typename T::value_type>);
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

  template <class T>
  static void container_swap(T& y, test::random_values<T> const& v)
  {
    T x(v.begin(), v.end());
    if (boost::allocator_propagate_on_container_swap<
          typename T::allocator_type>::type::value ||
        x.get_allocator() == y.get_allocator()) {
      y.swap(x);
    }
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
  }

#ifndef BOOST_UNORDERED_FOA_TESTS
  template <class T> static void buckets(T& y, test::random_values<T> const& v)
  {
    (void)y.begin(0);
    (void)y.end(0);
    (void)y.bucket_count();
    (void)y.max_bucket_count();
    (void)y.bucket_size(0);
    (void)y.bucket(get_key(*v.begin()));
  }
#endif

  template <class T>
  static void double_move_construct(T& y, test::random_values<T> const&)
  {
    T x = std::move(y);
    x.clear();
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
    BOOST_TEST_EQ(x.size(),
      static_cast<typename T::size_type>(std::distance(x.begin(), x.end())));
  }

  template <class T>
  static void double_move_assign(T& y, test::random_values<T> const&)
  {
    T x;
    x = std::move(y);
    x.clear();
    BOOST_TEST_EQ(y.size(),
      static_cast<typename T::size_type>(std::distance(y.begin(), y.end())));
    BOOST_TEST_EQ(x.size(),
      static_cast<typename T::size_type>(std::distance(x.begin(), x.end())));
  }

  template <class T>
  static void post_move_tests(T* ptr, test::random_generator const& generator)
  {
    // clang-format off
    void (*fps[])(T&, test::random_values<T> const&) = {
      insert_range<T>,
      insert_single<T>,
      insert_single_hint<T>,
      insert_or_assign<T>,
      insert_or_assign_hint<T>,
      try_emplace<T>,
      try_emplace_hint<T>,
      at<T>,
      index_operator<T>,
      clear<T>,
      capacity<T>,
      iterators<T>,
      erase_range<T>,
      erase_key<T>,
      lookup<T>,
      reserve<T>,
      copy_assignment<T>,
      move_assignment<T>,
      equal<T>,
      #ifndef BOOST_UNORDERED_FOA_TESTS
      extract<T>,
      buckets<T>,
      #endif
      merge<T>,
      erase_with_pred<T>,
      container_swap<T>,
      double_move_construct<T>,
      double_move_assign<T>
    };
    // clang-format on

    std::size_t const len = (sizeof(fps) / sizeof(*(fps)));

    for (std::size_t i = 0; i < len; ++i) {
      test::check_instances check_;

      test::random_values<T> const v(1000, generator);
      test::object_count count;
      T y(create(v, count));

      unsigned num_allocs = test::detail::tracker.count_allocations;
      (void)num_allocs;

      T x(std::move(y));

      BOOST_TEST(y.empty());
      BOOST_TEST(y.begin() == y.end());

#ifdef BOOST_UNORDERED_FOA_TESTS
      {
        using allocator_type = typename T::allocator_type;
        using value_type =
          typename boost::allocator_value_type<allocator_type>::type;
        using pointer = typename boost::allocator_pointer<allocator_type>::type;
        if (std::is_same<pointer, value_type*>::value) {
          BOOST_TEST_EQ(y.bucket_count(), 0u);
          BOOST_TEST_EQ(test::detail::tracker.count_allocations, num_allocs);
        }
      }
#else
      BOOST_TEST_EQ(y.bucket_count(), 0u);
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, num_allocs);
#endif

      fps[i](y, v);

      test::check_container(x, v);
      test::check_equivalent_keys(x);
    }

    for (std::size_t i = 0; i < len; ++i) {
      typename T::hasher hf(1);
      typename T::key_equal eq(1);
      typename T::allocator_type al1(1);
      typename T::allocator_type al2(2);

      test::check_instances check_;

      test::random_values<T> v(1000, generator);
      test::object_count count;
      T y(v.begin(), v.end(), 0, hf, eq, al1);
      T x(std::move(y), al2);

#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.empty());
      BOOST_TEST(y.begin() == y.end());
#else
      BOOST_TEST_NOT(y.empty());
      BOOST_TEST(y.begin() != y.end());
#endif

      fps[i](y, v);

      test::check_container(x, v);
      test::check_equivalent_keys(x);
    }

    for (std::size_t i = 0; i < len; ++i) {
      test::check_instances check_;

      test::random_values<T> v(1000, generator);
      test::object_count count;
      T y(create(v, count));

      unsigned num_allocs = test::detail::tracker.count_allocations;
      (void)num_allocs;

      T x(empty(ptr));
      x = std::move(y);

      BOOST_TEST(y.empty());
      BOOST_TEST(y.begin() == y.end());

#ifdef BOOST_UNORDERED_FOA_TESTS
      {
        using allocator_type = typename T::allocator_type;
        using value_type =
          typename boost::allocator_value_type<allocator_type>::type;
        using pointer = typename boost::allocator_pointer<allocator_type>::type;
        if (std::is_same<pointer, value_type*>::value) {
          BOOST_TEST_EQ(y.bucket_count(), 0u);
          BOOST_TEST_EQ(test::detail::tracker.count_allocations, num_allocs);
        }
      }
#else
      BOOST_TEST_EQ(y.bucket_count(), 0u);
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, num_allocs);
#endif

      fps[i](y, v);

      test::check_container(x, v);
      test::check_equivalent_keys(x);
    }

    for (std::size_t i = 0; i < len; ++i) {
      typename T::hasher hf(1);
      typename T::key_equal eq(1);
      typename T::allocator_type al1(1);
      typename T::allocator_type al2(2);

      test::check_instances check_;

      test::random_values<T> v(1000, generator);
      test::object_count count;
      T y(v.begin(), v.end(), 0, hf, eq, al1);

      unsigned num_allocs = test::detail::tracker.count_allocations;
      (void)num_allocs;

      T x(al2);
      x = std::move(y);

      bool b = boost::allocator_propagate_on_container_move_assignment<
        typename T::allocator_type>::type::value;
      if (b) {
        BOOST_TEST(y.empty());
        BOOST_TEST(y.begin() == y.end());
        BOOST_TEST_EQ(y.bucket_count(), 0u);
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, num_allocs);
      } else {
#ifdef BOOST_UNORDERED_FOA_TESTS
        BOOST_TEST(y.empty());
        BOOST_TEST(y.begin() == y.end());
#else
        BOOST_TEST_NOT(y.empty());
        BOOST_TEST(y.begin() != y.end());
#endif
      }

      fps[i](y, v);

      test::check_container(x, v);
      test::check_equivalent_keys(x);
    }
  }

  using test::default_generator;
  using test::generate_collisions;
  using test::limited_range;

#ifdef BOOST_UNORDERED_FOA_TESTS
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    std::allocator<std::pair<test::object const, test::object> > >*
    test_map_std_alloc;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    test::allocator2<test::object> >* test_set;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    test::allocator1<std::pair<test::object const, test::object> > >* test_map;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_move> >*
    test_set_prop_move;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::propagate_move> >* test_map_prop_move;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_move> >*
    test_set_no_prop_move;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::no_propagate_move> >* test_map_no_prop_move;

  boost::unordered_node_set<test::object, test::hash, test::equal_to,
    test::allocator2<test::object> >* test_node_set;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to,
    test::allocator1<std::pair<test::object const, test::object> > >*
    test_node_map;

  boost::unordered_node_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_move> >*
    test_node_set_prop_move;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::propagate_move> >* test_node_map_prop_move;

  boost::unordered_node_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_move> >*
    test_node_set_no_prop_move;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::no_propagate_move> >* test_node_map_no_prop_move;

  // clang-format off
  UNORDERED_TEST(post_move_tests,
    ((test_set)(test_map)(test_set_prop_move)(test_map_prop_move)
     (test_set_no_prop_move)(test_map_no_prop_move)
     (test_node_set)(test_node_map)
     (test_node_set_prop_move)(test_node_map_prop_move)
     (test_node_set_no_prop_move)(test_node_map_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
// clang-format on
#else
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    std::allocator<std::pair<test::object const, test::object> > >*
    test_map_std_alloc;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    test::allocator2<test::object> >* test_set;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_multiset;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    test::allocator1<std::pair<test::object const, test::object> > >* test_map;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to,
    test::allocator2<std::pair<test::object const, test::object> > >*
    test_multimap;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_move> >*
    test_set_prop_move;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_move> >*
    test_multiset_prop_move;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::propagate_move> >* test_map_prop_move;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::propagate_move> >* test_multimap_prop_move;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_move> >*
    test_set_no_prop_move;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_move> >*
    test_multiset_no_prop_move;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::no_propagate_move> >* test_map_no_prop_move;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<std::pair<test::object const, test::object>,
      test::no_propagate_move> >* test_multimap_no_prop_move;

  // clang-format off
  UNORDERED_TEST(post_move_tests,
    ((test_set)(test_multiset)(test_map)(test_multimap)(test_set_prop_move)(
      test_multiset_prop_move)(test_map_prop_move)(test_multimap_prop_move)(
      test_set_no_prop_move)(test_multiset_no_prop_move)(test_map_no_prop_move)(
      test_multimap_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
// clang-format on
#endif
} // namespace move_tests

RUN_TESTS()
