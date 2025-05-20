
// Copyright 2008-2009 Daniel James.
// Copyright 2022-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or move at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/equivalent.hpp"
#include "../helpers/invariants.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/test.hpp"
#include "../helpers/tracker.hpp"
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

  template <class T>
  void move_construct_tests1(T* ptr, test::random_generator const& generator)
  {
    typename T::hasher hf;
    typename T::key_equal eq;
    typename T::allocator_type al;

    {
      test::check_instances check_;

      T y(empty(ptr));
      BOOST_TEST(y.empty());
      BOOST_TEST(test::equivalent(y.hash_function(), hf));
      BOOST_TEST(test::equivalent(y.key_eq(), eq));
      BOOST_TEST(test::equivalent(y.get_allocator(), al));
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 1.0);
#endif
      test::check_equivalent_keys(y);

#ifdef BOOST_UNORDERED_FOA_TESTS
      using allocator_type = typename T::allocator_type;
      using value_type =
        typename boost::allocator_value_type<allocator_type>::type;
      using pointer = typename boost::allocator_pointer<allocator_type>::type;
      if (std::is_same<pointer, value_type*>::value) {
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }
#else
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
#endif
    }

    {
      test::check_instances check_;

      test::random_values<T> v(1000, generator);
      test::object_count count;
      T y(create(v, count));
#if defined(BOOST_HAS_NRVO)
      BOOST_TEST(count == test::global_object_count);
#endif
      test::check_container(y, v);
      test::check_equivalent_keys(y);
    }
  }

  template <class T>
  void move_assign_tests1(T* p, test::random_generator const& generator)
  {
    {
      test::check_instances check_;

      test::random_values<T> v(500, generator);
      test::object_count count;
      T y;
      y = create(v, count);
#if BOOST_UNORDERED_TEST_MOVING && defined(BOOST_HAS_NRVO)
      BOOST_TEST(count == test::global_object_count);
#endif
      test::check_container(y, v);
      test::check_equivalent_keys(y);
    }

    {
      test::random_values<T> v;

      T y;
      y = empty(p);

#ifdef BOOST_UNORDERED_FOA_TESTS
      using allocator_type = typename T::allocator_type;
      using value_type =
        typename boost::allocator_value_type<allocator_type>::type;
      using pointer = typename boost::allocator_pointer<allocator_type>::type;
      if (std::is_same<pointer, value_type*>::value) {
        BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
      }
#else
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
#endif
      test::check_container(y, v);
      test::check_equivalent_keys(y);
    }
  }

  template <class T>
  void move_construct_tests2(T*, test::random_generator const& generator)
  {
    typename T::hasher hf(1);
    typename T::key_equal eq(1);
    typename T::allocator_type al(1);
    typename T::allocator_type al2(2);

    test::object_count count;

    {
      test::check_instances check_;

      test::random_values<T> v(500, generator);
      T y(create(v, count, hf, eq, al, 0.5));
#if defined(BOOST_HAS_NRVO)
      BOOST_TEST(count == test::global_object_count);
#endif
      test::check_container(y, v);
      BOOST_TEST(test::equivalent(y.hash_function(), hf));
      BOOST_TEST(test::equivalent(y.key_eq(), eq));
      BOOST_TEST(test::equivalent(y.get_allocator(), al));
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 0.5); // Not necessarily required.
#endif
      test::check_equivalent_keys(y);
    }

    {
      test::check_instances check_;

      // TODO: To do this correctly requires the fancy new allocator
      // stuff.
      test::random_values<T> v(500, generator);
      T y(create(v, count, hf, eq, al, 2.0), al2);
      BOOST_TEST(count != test::global_object_count);
      test::check_container(y, v);
      BOOST_TEST(test::equivalent(y.hash_function(), hf));
      BOOST_TEST(test::equivalent(y.key_eq(), eq));
      BOOST_TEST(test::equivalent(y.get_allocator(), al2));
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 2.0); // Not necessarily required.
#endif
      test::check_equivalent_keys(y);
    }

    {
      test::check_instances check_;

      test::random_values<T> v(25, generator);
      T y(create(v, count, hf, eq, al, 1.0), al);
      BOOST_TEST(count == test::global_object_count);

      test::check_container(y, v);
      BOOST_TEST(test::equivalent(y.hash_function(), hf));
      BOOST_TEST(test::equivalent(y.key_eq(), eq));
      BOOST_TEST(test::equivalent(y.get_allocator(), al));
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 1.0); // Not necessarily required.
#endif
      test::check_equivalent_keys(y);
    }
  }

  template <class T>
  void move_assign_tests2(T*, test::random_generator const& generator)
  {
    typename T::hasher hf(1);
    typename T::key_equal eq(1);
    typename T::allocator_type al1(1);
    typename T::allocator_type al2(2);
    typedef typename T::allocator_type allocator_type;

    {
      test::random_values<T> v(500, generator);
      test::random_values<T> v2(0, generator);
      T y(v.begin(), v.end(), 0, hf, eq, al1);
      test::object_count count;
      y = create(v2, count, hf, eq, al2, 2.0);
      BOOST_TEST(y.empty());
      test::check_container(y, v2);
      test::check_equivalent_keys(y);
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 2.0);
#endif

#if defined(BOOST_HAS_NRVO)
      if (BOOST_UNORDERED_TEST_MOVING
            ? (bool)allocator_type::is_propagate_on_move
            : (bool)allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(y.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(y.get_allocator(), al1));
      }
#endif
    }

    {
      test::random_values<T> v(500, generator);
      test::object_count count;
      T y(0, hf, eq, al1);
      y = create(v, count, hf, eq, al2, 0.5);
#if defined(BOOST_HAS_NRVO)
      if (BOOST_UNORDERED_TEST_MOVING && allocator_type::is_propagate_on_move) {
        BOOST_TEST(count == test::global_object_count);
      }
#endif
      test::check_container(y, v);
      test::check_equivalent_keys(y);
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 0.5);
#endif

#if defined(BOOST_HAS_NRVO)
      if (BOOST_UNORDERED_TEST_MOVING
            ? (bool)allocator_type::is_propagate_on_move
            : (bool)allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(y.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(y.get_allocator(), al1));
      }
#endif
    }

    {
      test::random_values<T> v;
      T y(0, hf, eq, al1);
      T x(0, hf, eq, al2);
      x.max_load_factor(0.25);

#ifdef BOOST_UNORDERED_FOA_TESTS
      {
        using value_type =
          typename boost::allocator_value_type<allocator_type>::type;
        using pointer = typename boost::allocator_pointer<allocator_type>::type;
        if (std::is_same<pointer, value_type*>::value) {
          BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
        }
      }
#else
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
#endif

      y = std::move(x);

#ifdef BOOST_UNORDERED_FOA_TESTS
      {
        using value_type =
          typename boost::allocator_value_type<allocator_type>::type;
        using pointer = typename boost::allocator_pointer<allocator_type>::type;
        if (std::is_same<pointer, value_type*>::value) {
          BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
        }
      }
#else
      BOOST_TEST_EQ(test::detail::tracker.count_allocations, 0u);
#endif
      test::check_container(y, v);
      test::check_equivalent_keys(y);
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 0.25);
#endif

      if (BOOST_UNORDERED_TEST_MOVING
            ? (bool)allocator_type::is_propagate_on_move
            : (bool)allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(y.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(y.get_allocator(), al1));
      }
    }

    {
      test::check_instances check_;

      test::random_values<T> v(500, generator);
      T y(0, hf, eq, al1);

      T x(0, hf, eq, al2);
      x.max_load_factor(0.25);
      x.insert(v.begin(), v.end());

      test::object_count count = test::global_object_count;
      y = std::move(x);
      if (BOOST_UNORDERED_TEST_MOVING && allocator_type::is_propagate_on_move) {
        BOOST_TEST(count == test::global_object_count);
      }
      test::check_container(y, v);
      test::check_equivalent_keys(y);
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 0.25);
#endif

      if (BOOST_UNORDERED_TEST_MOVING
            ? (bool)allocator_type::is_propagate_on_move
            : (bool)allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(y.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(y.get_allocator(), al1));
      }
    }

    {
      test::check_instances check_;

      test::random_values<T> v1(1000, generator);
      test::random_values<T> v2(200, generator);

      T x(0, hf, eq, al2);
      x.max_load_factor(0.5);
      x.insert(v2.begin(), v2.end());

      test::object_count count1 = test::global_object_count;

      T y(v1.begin(), v1.end(), 0, hf, eq, al1);
      y = std::move(x);

      test::object_count count2 = test::global_object_count;

      if (BOOST_UNORDERED_TEST_MOVING && allocator_type::is_propagate_on_move) {
        BOOST_TEST(count1.instances == test::global_object_count.instances);
        BOOST_TEST(
          count2.constructions == test::global_object_count.constructions);
      }

      test::check_container(y, v2);
      test::check_equivalent_keys(y);
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST(y.max_load_factor() == 0.875);
#else
      BOOST_TEST(y.max_load_factor() == 0.5);
#endif

      if (BOOST_UNORDERED_TEST_MOVING
            ? (bool)allocator_type::is_propagate_on_move
            : (bool)allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(y.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(y.get_allocator(), al1));
      }
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

  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to,
    std::allocator<std::pair<test::object const, test::object> > >*
    test_node_map_std_alloc;

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
  UNORDERED_TEST(move_construct_tests1,
    ((test_map_std_alloc)(test_set)(test_map)
     (test_set_prop_move)(test_map_prop_move)
     (test_set_no_prop_move)(test_map_no_prop_move)
     (test_node_map_std_alloc)(test_node_set)(test_node_map)
     (test_node_set_prop_move)(test_node_map_prop_move)
     (test_node_set_no_prop_move)(test_node_map_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(move_assign_tests1,
    ((test_map_std_alloc)(test_set)(test_map)
     (test_set_prop_move)(test_map_prop_move)
     (test_set_no_prop_move)(test_map_no_prop_move)
     (test_node_map_std_alloc)(test_node_set)(test_node_map)
     (test_node_set_prop_move)(test_node_map_prop_move)
     (test_node_set_no_prop_move)(test_node_map_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(move_construct_tests2,
    ((test_set)(test_map)
     (test_set_prop_move)(test_map_prop_move)
     (test_set_no_prop_move)(test_map_no_prop_move)
     (test_node_set)(test_node_map)
     (test_node_set_prop_move)(test_node_map_prop_move)
     (test_node_set_no_prop_move)(test_node_map_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(move_assign_tests2,
    ((test_set)(test_map)
     (test_set_prop_move)(test_map_prop_move)
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

  UNORDERED_TEST(move_construct_tests1,
    ((test_map_std_alloc)(test_set)(test_multiset)(test_map)(test_multimap)(test_set_prop_move)(test_multiset_prop_move)(test_map_prop_move)(test_multimap_prop_move)(test_set_no_prop_move)(test_multiset_no_prop_move)(test_map_no_prop_move)(test_multimap_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(move_assign_tests1,
    ((test_map_std_alloc)(test_set)(test_multiset)(test_map)(test_multimap)(test_set_prop_move)(test_multiset_prop_move)(test_map_prop_move)(test_multimap_prop_move)(test_set_no_prop_move)(test_multiset_no_prop_move)(test_map_no_prop_move)(test_multimap_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(move_construct_tests2,
    ((test_set)(test_multiset)(test_map)(test_multimap)(test_set_prop_move)(test_multiset_prop_move)(test_map_prop_move)(test_multimap_prop_move)(test_set_no_prop_move)(test_multiset_no_prop_move)(test_map_no_prop_move)(test_multimap_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
  UNORDERED_TEST(move_assign_tests2,
    ((test_set)(test_multiset)(test_map)(test_multimap)(test_set_prop_move)(test_multiset_prop_move)(test_map_prop_move)(test_multimap_prop_move)(test_set_no_prop_move)(test_multiset_no_prop_move)(test_map_no_prop_move)(test_multimap_no_prop_move))(
      (default_generator)(generate_collisions)(limited_range)))
#endif
} // namespace move_tests

RUN_TESTS()
