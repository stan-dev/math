
// Copyright 2006-2009 Daniel James.
// Copyright 2022-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"
#include "../objects/test.hpp"
#include "../objects/cxx11_allocator.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/tracker.hpp"
#include "../helpers/equivalent.hpp"

#if defined(BOOST_MSVC)
#pragma warning(disable : 4127) // conditional expression is constant
#endif

#if defined(__clang__) && defined(__has_warning)
#if __has_warning("-Wself-assign-overloaded")
#pragma clang diagnostic ignored "-Wself-assign-overloaded"
#endif
#endif

namespace assign_tests {

  test::seed_t initialize_seed(96785);

  template <class T> void assign_tests1(T*, test::random_generator generator)
  {
    typename T::hasher hf;
    typename T::key_equal eq;

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests1.1\n";
    {
      test::check_instances check_;

      T x;
      x = x;
      BOOST_TEST(x.empty());
      BOOST_TEST(test::equivalent(x.hash_function(), hf));
      BOOST_TEST(test::equivalent(x.key_eq(), eq));
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests1.2\n";
    {
      test::check_instances check_;

      test::random_values<T> v(1000, generator);
      T x(v.begin(), v.end());

      test::ordered<T> tracker = test::create_ordered(x);
      tracker.insert_range(v.begin(), v.end());

      x = x;
      tracker.compare(x);

      T y;
      y.max_load_factor(x.max_load_factor() / 20);
      float mlf = x.max_load_factor();
      y = x;
      tracker.compare(x);
      tracker.compare(y);
      BOOST_TEST(x.max_load_factor() == mlf);
      BOOST_TEST(y.max_load_factor() == mlf);
      BOOST_TEST(y.load_factor() <= y.max_load_factor());
    }
  }

  template <class T> void assign_tests2(T*, test::random_generator generator)
  {
    typename T::hasher hf1(1);
    typename T::hasher hf2(2);
    typename T::key_equal eq1(1);
    typename T::key_equal eq2(2);
    typename T::allocator_type al1(1);
    typename T::allocator_type al2(2);

    typedef typename T::allocator_type allocator_type;

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests2.0 - empty container\n";
    {
      test::check_instances check_;

      T x1(0, hf1, eq1);
      T x2(0, hf2, eq2);
      x2 = x1;
      BOOST_TEST(test::equivalent(x1.hash_function(), hf1));
      BOOST_TEST(test::equivalent(x1.key_eq(), eq1));
      BOOST_TEST(test::equivalent(x2.hash_function(), hf1));
      BOOST_TEST(test::equivalent(x2.key_eq(), eq1));
      test::check_container(x1, x2);
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests2.1\n";
    {
      test::check_instances check_;

      test::random_values<T> v(1000, generator);
      T x1(v.begin(), v.end(), 0, hf1, eq1);
      T x2(0, hf2, eq2);
      x2 = x1;
      BOOST_TEST(test::equivalent(x1.hash_function(), hf1));
      BOOST_TEST(test::equivalent(x1.key_eq(), eq1));
      BOOST_TEST(test::equivalent(x2.hash_function(), hf1));
      BOOST_TEST(test::equivalent(x2.key_eq(), eq1));
      test::check_container(x1, v);
      test::check_container(x2, v);
      BOOST_TEST(x2.load_factor() <= x2.max_load_factor());
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests2.1a\n";
    {
      test::check_instances check_;

      test::random_values<T> v1(0, generator);
      test::random_values<T> v2(1000, generator);
      T x1(0, hf2, eq2);
      T x2(v2.begin(), v2.end(), 0, hf1, eq1);
      x2 = x1;
      BOOST_TEST(test::equivalent(x1.hash_function(), hf2));
      BOOST_TEST(test::equivalent(x1.key_eq(), eq2));
      BOOST_TEST(test::equivalent(x2.hash_function(), hf2));
      BOOST_TEST(test::equivalent(x2.key_eq(), eq2));
      test::check_container(x1, v1);
      test::check_container(x2, v1);
      BOOST_TEST(x2.load_factor() <= x2.max_load_factor());
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests2.2\n";
    {
      test::check_instances check_;

      test::random_values<T> v1(100, generator), v2(100, generator);
      T x1(v1.begin(), v1.end(), 0, hf1, eq1, al1);
      T x2(v2.begin(), v2.end(), 0, hf2, eq2, al2);
      x2 = x1;
      BOOST_TEST(test::equivalent(x2.hash_function(), hf1));
      BOOST_TEST(test::equivalent(x2.key_eq(), eq1));
      if (allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(x2.get_allocator(), al1));
        BOOST_TEST(!test::equivalent(x2.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(x2.get_allocator(), al2));
        BOOST_TEST(!test::equivalent(x2.get_allocator(), al1));
      }
      test::check_container(x1, v1);
      test::check_container(x2, v1);
      BOOST_TEST(x2.load_factor() <= x2.max_load_factor());
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests2.3\n";
    {
      test::check_instances check_;

      test::random_values<T> v1(100, generator), v2(1000, generator);
      T x1(v1.begin(), v1.end(), 0, hf1, eq1, al1);
      T x2(v2.begin(), v2.end(), 0, hf2, eq2, al2);
      x2 = x1;
      BOOST_TEST(test::equivalent(x2.hash_function(), hf1));
      BOOST_TEST(test::equivalent(x2.key_eq(), eq1));
      if (allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(x2.get_allocator(), al1));
        BOOST_TEST(!test::equivalent(x2.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(x2.get_allocator(), al2));
        BOOST_TEST(!test::equivalent(x2.get_allocator(), al1));
      }
      test::check_container(x1, v1);
      test::check_container(x2, v1);
      BOOST_TEST(x2.load_factor() <= x2.max_load_factor());
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "assign_tests2.4\n";
    {
      test::check_instances check_;

      test::random_values<T> v1(1000, generator), v2(100, generator);
      T x1(v1.begin(), v1.end(), 0, hf1, eq1, al1);
      T x2(v2.begin(), v2.end(), 0, hf2, eq2, al2);
      x2 = x1;
      BOOST_TEST(test::equivalent(x2.hash_function(), hf1));
      BOOST_TEST(test::equivalent(x2.key_eq(), eq1));
      if (allocator_type::is_propagate_on_assign) {
        BOOST_TEST(test::equivalent(x2.get_allocator(), al1));
        BOOST_TEST(!test::equivalent(x2.get_allocator(), al2));
      } else {
        BOOST_TEST(test::equivalent(x2.get_allocator(), al2));
        BOOST_TEST(!test::equivalent(x2.get_allocator(), al1));
      }
      test::check_container(x1, v1);
      test::check_container(x2, v1);
      BOOST_TEST(x2.load_factor() <= x2.max_load_factor());
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "gh205\n";
    {
      // https://github.com/boostorg/unordered/issues/205
      {
        // A=B
        // assign no-allocated to empty-but-allocated
        test::check_instances check_;

        test::random_values<T> v1(1, generator);
        T x1(v1.begin(), v1.end(), 0, hf1, eq1, al1);
        T x2(0, hf2, eq2, al2);
        x1.clear();
        x1 = x2;
        BOOST_TEST(x1.empty());
        BOOST_TEST(x1.begin() == x1.end());
        BOOST_TEST(x2.empty());
        BOOST_TEST(x2.begin() == x2.end());
      }

      {
        // B=A
        // assign empty-but-allocated to no-allocated
        test::check_instances check_;

        test::random_values<T> v1(1, generator);
        T x1(v1.begin(), v1.end(), 0, hf1, eq1, al1);
        T x2(0, hf2, eq2, al2);
        x1.clear();
        x2 = x1;
        BOOST_TEST(x2.empty());
        BOOST_TEST(x2.begin() == x2.end());
        BOOST_TEST(x1.empty());
        BOOST_TEST(x1.begin() == x1.end());
      }

      {
        // A=A
        // assign empty-but-allocated to empty-but-allocated
        test::check_instances check_;

        test::random_values<T> v1(1, generator);
        test::random_values<T> v2(1, generator);
        T x1(v1.begin(), v1.end(), 0, hf1, eq1, al1);
        T x2(v2.begin(), v2.end(), 0, hf2, eq2, al2);
        x1.clear();
        x2.clear();
        x1 = x2;
        BOOST_TEST(x1.empty());
        BOOST_TEST(x1.begin() == x1.end());
        BOOST_TEST(x2.empty());
        BOOST_TEST(x2.begin() == x2.end());
      }

      {
        // B=B
        // assign no-allocated to no-allocated
        test::check_instances check_;

        T x1(0, hf1, eq1, al1);
        T x2(0, hf2, eq2, al2);
        x1 = x2;
        BOOST_TEST(x2.empty());
        BOOST_TEST(x2.begin() == x2.end());
        BOOST_TEST(x1.empty());
        BOOST_TEST(x1.begin() == x1.end());
      }

#ifdef BOOST_UNORDERED_FOA_TESTS
      {
        // check that optimized copying preserves the anti-drift mechanism
        test::check_instances check_;

        test::random_values<T> v1(2000, generator);
        test::random_values<T> v2(2000, generator);
        T x1(v1.begin(), v1.end(), 0, hf1, eq1, al1);
        T x2(v2.begin(), v2.end(), 0, hf2, eq2, al2);

        auto ml = x1.max_load();

        while (ml == x1.max_load() && x1.size()) {
          x1.erase(x1.begin());
        };

        BOOST_TEST_EQ(x1.bucket_count(), x2.bucket_count());

        x2 = x1;

        if (x1.max_load() != x2.max_load()) {
          BOOST_TEST(boost::allocator_propagate_on_container_copy_assignment<
            typename T::allocator_type>::type::value);
        }
      }
#endif
    }
  }

  using test::default_generator;
  using test::generate_collisions;
  using test::limited_range;

  template <typename T> bool is_propagate(T*)
  {
    return T::allocator_type::is_propagate_on_assign;
  }

#ifdef BOOST_UNORDERED_FOA_TESTS
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to, std::allocator<test::object> >* test_map_std_alloc;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to, std::allocator<test::object> >* test_node_map_std_alloc;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_set;
  boost::unordered_node_set<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_node_set;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to, test::allocator2<test::object> >* test_map;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to, test::allocator2<test::object> >* test_node_map;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_set_prop_assign;
  boost::unordered_node_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_node_set_prop_assign;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_map_prop_assign;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_node_map_prop_assign;

  boost::unordered_flat_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_set_no_prop_assign;
  boost::unordered_node_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_node_set_no_prop_assign;
  boost::unordered_flat_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_map_no_prop_assign;
  boost::unordered_node_map<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_node_map_no_prop_assign;

  UNORDERED_AUTO_TEST (check_traits) {
    BOOST_TEST(!is_propagate(test_set));
    BOOST_TEST(is_propagate(test_set_prop_assign));
    BOOST_TEST(!is_propagate(test_set_no_prop_assign));

    BOOST_TEST(!is_propagate(test_node_set));
    BOOST_TEST(is_propagate(test_node_set_prop_assign));
    BOOST_TEST(!is_propagate(test_node_set_no_prop_assign));
  }

  UNORDERED_TEST(assign_tests1,
    ((test_map_std_alloc)(test_node_map_std_alloc)
     (test_set)(test_node_set)
     (test_map)(test_node_map)
     (test_set_prop_assign)(test_node_set_prop_assign)
     (test_map_prop_assign)(test_node_map_prop_assign)
     (test_set_no_prop_assign)(test_node_set_no_prop_assign)
     (test_map_no_prop_assign)(test_node_map_no_prop_assign))(
      (default_generator)(generate_collisions)(limited_range)))

  UNORDERED_TEST(assign_tests2,
    ((test_set)(test_node_set)
     (test_map)(test_node_map)
     (test_set_prop_assign)(test_node_set_prop_assign)
     (test_map_prop_assign)(test_node_map_prop_assign)
     (test_set_no_prop_assign)(test_node_set_no_prop_assign)
     (test_map_no_prop_assign)(test_node_map_no_prop_assign))(
      (default_generator)(generate_collisions)(limited_range)))
#else
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    std::allocator<test::object> >* test_map_std_alloc;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_set;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::allocator2<test::object> >* test_multiset;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    test::allocator2<test::object> >* test_map;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to, test::allocator1<test::object> >* test_multimap;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_set_prop_assign;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_multiset_prop_assign;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_map_prop_assign;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<test::object, test::propagate_assign> >*
    test_multimap_prop_assign;

  boost::unordered_set<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_set_no_prop_assign;
  boost::unordered_multiset<test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_multiset_no_prop_assign;
  boost::unordered_map<test::object, test::object, test::hash, test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_map_no_prop_assign;
  boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to,
    test::cxx11_allocator<test::object, test::no_propagate_assign> >*
    test_multimap_no_prop_assign;

  UNORDERED_AUTO_TEST (check_traits) {
    BOOST_TEST(!is_propagate(test_set));
    BOOST_TEST(is_propagate(test_set_prop_assign));
    BOOST_TEST(!is_propagate(test_set_no_prop_assign));
  }

  UNORDERED_TEST(assign_tests1,
    ((test_map_std_alloc)(test_set)(test_multiset)(test_map)(test_multimap)(
      test_set_prop_assign)(test_multiset_prop_assign)(test_map_prop_assign)(
      test_multimap_prop_assign)(test_set_no_prop_assign)(
      test_multiset_no_prop_assign)(test_map_no_prop_assign)(
      test_multimap_no_prop_assign))(
      (default_generator)(generate_collisions)(limited_range)))

  UNORDERED_TEST(
    assign_tests2, ((test_set)(test_multiset)(test_map)(test_multimap)(
                     test_set_prop_assign)(test_multiset_prop_assign)(
                     test_map_prop_assign)(test_multimap_prop_assign)(
                     test_set_no_prop_assign)(test_multiset_no_prop_assign)(
                     test_map_no_prop_assign)(test_multimap_no_prop_assign))(
      (default_generator)(generate_collisions)(limited_range)))
#endif

  UNORDERED_AUTO_TEST (assign_default_initializer_list) {
    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Initializer List Tests\n";
    std::initializer_list<std::pair<int const, int> > init;
#ifdef BOOST_UNORDERED_FOA_TESTS
    boost::unordered_flat_map<int, int> x1;
    boost::unordered_node_map<int, int> x2;
    x2[25] = 3;
    x2[16] = 10;
    BOOST_TEST(!x2.empty());
    x2 = init;
    BOOST_TEST(x2.empty());
#else
    boost::unordered_map<int, int> x1;
#endif
    x1[25] = 3;
    x1[16] = 10;
    BOOST_TEST(!x1.empty());
    x1 = init;
    BOOST_TEST(x1.empty());
  }

  UNORDERED_AUTO_TEST (assign_initializer_list) {
    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Initializer List Tests\n";

#ifdef BOOST_UNORDERED_FOA_TESTS
    boost::unordered_flat_set<int> x;
    boost::unordered_node_set<int> y;
    y.insert(10);
    y.insert(20);
    y = {1, 2, -10};
    BOOST_TEST(y.find(10) == y.end());
    BOOST_TEST(y.find(-10) != y.end());
#else
    boost::unordered_set<int> x;
#endif
    x.insert(10);
    x.insert(20);
    x = {1, 2, -10};
    BOOST_TEST(x.find(10) == x.end());
    BOOST_TEST(x.find(-10) != x.end());
  }
}

RUN_TESTS()
