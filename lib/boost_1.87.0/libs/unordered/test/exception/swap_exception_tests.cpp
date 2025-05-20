
// Copyright 2006-2009 Daniel James.
// Copyright 2022 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_ENABLE_ASSERT_HANDLER
#include <boost/assert.hpp>

#include "./containers.hpp"

#include "../helpers/invariants.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/tracker.hpp"
#include "../objects/test.hpp"

#include <sstream>

namespace boost {
  void assertion_failed(
    char const* expr, char const* function, char const* file, long line)
  {
    std::stringstream ss;
    ss << expr << "\nin " << function << " failed at : " << file << ", line "
       << line;

    throw std::runtime_error(ss.str());
  }

  void assertion_failed_msg(char const* expr, char const* msg,
    char const* function, char const* file, long line)
  {
    std::stringstream ss;
    ss << expr << "\nin " << function << " failed at : " << file << ", line "
       << line << "\n"
       << msg;

    throw std::runtime_error(ss.str());
  }
} // namespace boost

#if defined(BOOST_MSVC)
#pragma warning(disable : 4512) // assignment operator could not be generated
#endif

test::seed_t initialize_seed(9387);

template <class T> struct self_swap_base : public test::exception_base
{
  test::random_values<T> values;
  self_swap_base(std::size_t count = 0) : values(count, test::limited_range) {}

  typedef T data_type;
  T init() const { return T(values.begin(), values.end()); }

  void run(T& x) const
  {
    x.swap(x);

    DISABLE_EXCEPTIONS;
    test::check_container(x, this->values);
    test::check_equivalent_keys(x);
  }

  void check BOOST_PREVENT_MACRO_SUBSTITUTION(T const& x) const
  {
    (void)x;

    BOOST_ERROR("An exception leaked when it should not have. Allocator "
                "equality assertion must precede all other ops");
  }
};

template <class T> struct self_swap_test1 : self_swap_base<T>
{
};

template <class T> struct self_swap_test2 : self_swap_base<T>
{
  self_swap_test2() : self_swap_base<T>(100) {}
};

template <class T> struct swap_base : public test::exception_base
{
  const test::random_values<T> x_values, y_values;
  const T initial_x, initial_y;

  typedef typename T::hasher hasher;
  typedef typename T::key_equal key_equal;
  typedef typename T::allocator_type allocator_type;

  swap_base(unsigned int count1, unsigned int count2, int tag1, int tag2)
      : x_values(count1, test::limited_range),
        y_values(count2, test::limited_range),
        initial_x(x_values.begin(), x_values.end(), 0, hasher(tag1),
          key_equal(tag1), allocator_type(tag1)),
        initial_y(y_values.begin(), y_values.end(), 0, hasher(tag2),
          key_equal(tag2),
          allocator_type(T::allocator_type::propagate_on_container_swap::value
                           ? tag2
                           : tag1))
  {
  }

  struct data_type
  {
    data_type(T const& x_, T const& y_) : x(x_), y(y_) {}

    T x, y;
  };

  data_type init() const { return data_type(initial_x, initial_y); }

  void run(data_type& d) const
  {
    try {
      d.x.swap(d.y);
    } catch (std::runtime_error&) {
    }

    DISABLE_EXCEPTIONS;
    test::check_container(d.x, this->y_values);
    test::check_equivalent_keys(d.x);
    test::check_container(d.y, this->x_values);
    test::check_equivalent_keys(d.y);
  }

  void check BOOST_PREVENT_MACRO_SUBSTITUTION(data_type const& d) const
  {
    std::string scope(test::scope);

    // TODO: In C++11 exceptions are only allowed in the swap function.
    BOOST_TEST(scope == "hash::hash(hash)" ||
               scope == "hash::operator=(hash)" ||
               scope == "equal_to::equal_to(equal_to)" ||
               scope == "equal_to::operator=(equal_to)");

    test::check_equivalent_keys(d.x);
    test::check_equivalent_keys(d.y);
  }
};

template <class T> struct swap_test1 : swap_base<T>
{
  swap_test1() : swap_base<T>(0, 0, 0, 0) {}
};

template <class T> struct swap_test2 : swap_base<T>
{
  swap_test2() : swap_base<T>(60, 0, 0, 0) {}
};

template <class T> struct swap_test3 : swap_base<T>
{
  swap_test3() : swap_base<T>(0, 60, 0, 0) {}
};

template <class T> struct swap_test4 : swap_base<T>
{
  swap_test4() : swap_base<T>(10, 10, 1, 2) {}
};

template <class T> struct unequal_alloc_swap_base : public test::exception_base
{
  const test::random_values<T> x_values, y_values;
  const T initial_x, initial_y;

  typedef typename T::hasher hasher;
  typedef typename T::key_equal key_equal;
  typedef typename T::allocator_type allocator_type;

  unequal_alloc_swap_base(unsigned int count1, unsigned int count2)
      : x_values(count1, test::limited_range),
        y_values(count2, test::limited_range),
        initial_x(x_values.begin(), x_values.end(), 0, allocator_type(1337)),
        initial_y(y_values.begin(), y_values.end(), 0, allocator_type(7331))
  {
  }

  struct data_type
  {
    data_type(T const& x_, T const& y_) : x(x_), y(y_) {}

    T x, y;
  };

  data_type init() const { return data_type(initial_x, initial_y); }

  void run(data_type& d) const
  {
    bool assert_threw = false;

    BOOST_TEST(d.x.get_allocator() != d.y.get_allocator());

    try {
      d.x.swap(d.y);
    } catch (std::runtime_error&) {
      assert_threw = true;
    }

    DISABLE_EXCEPTIONS;
    BOOST_TEST(assert_threw);
    test::check_container(d.x, this->x_values);
    test::check_equivalent_keys(d.x);
    test::check_container(d.y, this->y_values);
    test::check_equivalent_keys(d.y);
  }

  void check BOOST_PREVENT_MACRO_SUBSTITUTION(data_type const& d) const
  {
    std::string scope(test::scope);

    // TODO: In C++11 exceptions are only allowed in the swap function.
    BOOST_TEST(scope == "hash::hash(hash)" ||
               scope == "hash::operator=(hash)" ||
               scope == "equal_to::equal_to(equal_to)" ||
               scope == "equal_to::operator=(equal_to)");

    test::check_equivalent_keys(d.x);
    test::check_equivalent_keys(d.y);
  }
};

template <class T> struct unequal_alloc_swap_test1 : unequal_alloc_swap_base<T>
{
  unequal_alloc_swap_test1() : unequal_alloc_swap_base<T>(0, 0) {}
};

template <class T> struct unequal_alloc_swap_test2 : unequal_alloc_swap_base<T>
{
  unequal_alloc_swap_test2() : unequal_alloc_swap_base<T>(0, 10) {}
};

template <class T> struct unequal_alloc_swap_test3 : unequal_alloc_swap_base<T>
{
  unequal_alloc_swap_test3() : unequal_alloc_swap_base<T>(10, 0) {}
};

template <class T> struct unequal_alloc_swap_test4 : unequal_alloc_swap_base<T>
{
  unequal_alloc_swap_test4() : unequal_alloc_swap_base<T>(10, 10) {}
};

#if defined(BOOST_UNORDERED_FOA_TESTS)

using unordered_flat_set = boost::unordered_flat_set<int, boost::hash<int>,
  std::equal_to<int>, test::allocator1<int> >;
using unordered_flat_map = boost::unordered_flat_map<int, int, boost::hash<int>,
  std::equal_to<int>, test::allocator1<std::pair<int const, int> > >;
using unordered_node_set = boost::unordered_node_set<int, boost::hash<int>,
  std::equal_to<int>, test::allocator1<int> >;
using unordered_node_map = boost::unordered_node_map<int, int, boost::hash<int>,
  std::equal_to<int>, test::allocator1<std::pair<int const, int> > >;

#define SWAP_CONTAINER_SEQ \
 (unordered_flat_set)(unordered_flat_map) \
 (unordered_node_set)(unordered_node_map)

#else

typedef boost::unordered_set<int, boost::hash<int>, std::equal_to<int>,
  test::allocator1<int> >
  unordered_set;
typedef boost::unordered_map<int, int, boost::hash<int>, std::equal_to<int>,
  test::allocator1<std::pair<int const, int> > >
  unordered_map;
typedef boost::unordered_multiset<int, boost::hash<int>, std::equal_to<int>,
  test::allocator1<int> >
  unordered_multiset;
typedef boost::unordered_multimap<int, int, boost::hash<int>,
  std::equal_to<int>, test::allocator1<std::pair<int const, int> > >
  unordered_multimap;

#define SWAP_CONTAINER_SEQ                                                     \
  (unordered_set)(unordered_map)(unordered_multiset)(unordered_multimap)
#endif

// FOA containers deliberately choose to not offer the strong exception
// guarantee so we can't reliably test what happens if swapping one of the data
// members throws
//
// clang-format off
#if !defined(BOOST_UNORDERED_FOA_TESTS)
EXCEPTION_TESTS(
  (self_swap_test1)(self_swap_test2)
  (swap_test1)(swap_test2)(swap_test3)(swap_test4),
  CONTAINER_SEQ)
#endif

// want to prove that when assertions are defined as throwing operations that we
// uphold invariants
EXCEPTION_TESTS(
  (unequal_alloc_swap_test1)(unequal_alloc_swap_test2)
  (unequal_alloc_swap_test3)(unequal_alloc_swap_test4),
  SWAP_CONTAINER_SEQ)
// clang-format on

RUN_TESTS()
