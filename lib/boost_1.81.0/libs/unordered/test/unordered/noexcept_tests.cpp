
// Copyright 2013 Daniel James.
// Copyright 2022 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"
#include "../helpers/fwd.hpp"

#if defined(BOOST_MSVC)
#pragma warning(push)
// conditional expression is constant
#pragma warning(disable : 4127)
#endif

namespace noexcept_tests {
  // Test the noexcept is set correctly for the move constructor.

  struct hash_possible_exception : boost::hash<int>
  {
    hash_possible_exception(hash_possible_exception const&) {}
    hash_possible_exception& operator=(hash_possible_exception const&)
    {
      return *this;
    }
  };

  struct equal_to_possible_exception : std::equal_to<int>
  {
    equal_to_possible_exception(equal_to_possible_exception const&) {}
    equal_to_possible_exception& operator=(equal_to_possible_exception const&)
    {
      return *this;
    }
  };

  // Test that the move constructor does actually move without throwing
  // an exception when it claims to.

  struct test_exception
  {
  };

  bool throwing_test_exception = false;
  void test_throw(char const* name)
  {
    if (throwing_test_exception) {
      BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Throw exception in: " << name
                                     << std::endl;
      throw test_exception();
    }
  }

  template <bool nothrow_move_construct, bool nothrow_move_assign,
    bool nothrow_swap>
  class hash_nothrow : boost::hash<int>
  {
    BOOST_COPYABLE_AND_MOVABLE(hash_nothrow)

    typedef boost::hash<int> base;

  public:
    hash_nothrow(BOOST_RV_REF(hash_nothrow))
      BOOST_NOEXCEPT_IF(nothrow_move_construct)
    {
      if (!nothrow_move_construct) {
        test_throw("Move Constructor");
      }
    }

    hash_nothrow() { test_throw("Constructor"); }
    hash_nothrow(hash_nothrow const&) { test_throw("Copy"); }
    hash_nothrow& operator=(BOOST_COPY_ASSIGN_REF(hash_nothrow))
    {
      test_throw("Assign");
      return *this;
    }
    hash_nothrow& operator=(BOOST_RV_REF(hash_nothrow))
      BOOST_NOEXCEPT_IF(nothrow_move_assign)
    {
      if (!nothrow_move_assign) {
        test_throw("Move Assign");
      }
      return *this;
    }
    std::size_t operator()(int x) const
    {
      test_throw("Operator");
      return static_cast<base const&>(*this)(x);
    }
    friend void swap(hash_nothrow&, hash_nothrow&)
      BOOST_NOEXCEPT_IF(nothrow_swap)
    {
      if (!nothrow_swap) {
        test_throw("Swap");
      }
    }
  };

  typedef hash_nothrow<true, false, false> hash_nothrow_move_construct;
  typedef hash_nothrow<false, true, false> hash_nothrow_move_assign;
  typedef hash_nothrow<false, false, true> hash_nothrow_swap;

  template <bool nothrow_move_construct, bool nothrow_move_assign,
    bool nothrow_swap>
  class equal_to_nothrow
  {
    BOOST_COPYABLE_AND_MOVABLE(equal_to_nothrow)

    typedef boost::hash<int> base;

  public:
    equal_to_nothrow(BOOST_RV_REF(equal_to_nothrow))
      BOOST_NOEXCEPT_IF(nothrow_move_construct)
    {
      if (!nothrow_move_construct) {
        test_throw("Move Constructor");
      }
    }

    equal_to_nothrow() { test_throw("Constructor"); }
    equal_to_nothrow(equal_to_nothrow const&) { test_throw("Copy"); }
    equal_to_nothrow& operator=(BOOST_COPY_ASSIGN_REF(equal_to_nothrow))
    {
      test_throw("Assign");
      return *this;
    }
    equal_to_nothrow& operator=(BOOST_RV_REF(equal_to_nothrow))
      BOOST_NOEXCEPT_IF(nothrow_move_assign)
    {
      if (!nothrow_move_assign) {
        test_throw("Move Assign");
      }
      return *this;
    }
    std::size_t operator()(int x, int y) const
    {
      test_throw("Operator");
      return x == y;
    }
    friend void swap(equal_to_nothrow&, equal_to_nothrow&)
      BOOST_NOEXCEPT_IF(nothrow_swap)
    {
      if (!nothrow_swap) {
        test_throw("Swap");
      }
    }
  };

  typedef equal_to_nothrow<true, false, false> equal_to_nothrow_move_construct;
  typedef equal_to_nothrow<false, true, false> equal_to_nothrow_move_assign;
  typedef equal_to_nothrow<false, false, true> equal_to_nothrow_swap;

  bool have_is_nothrow_move = false;
  bool have_is_nothrow_move_assign = false;
  bool have_is_nothrow_swap = false;

  UNORDERED_AUTO_TEST (check_is_nothrow_move) {
    BOOST_TEST(
      !boost::is_nothrow_move_constructible<hash_possible_exception>::value);
    BOOST_TEST(
      !boost::is_nothrow_move_assignable<hash_possible_exception>::value);
    BOOST_TEST(!boost::is_nothrow_swappable<hash_possible_exception>::value);
    BOOST_TEST((!boost::is_nothrow_move_constructible<
                equal_to_nothrow<false, false, false> >::value));
    BOOST_TEST((!boost::is_nothrow_move_assignable<
                equal_to_nothrow<false, false, false> >::value));
    BOOST_TEST((!boost::is_nothrow_swappable<
                equal_to_nothrow<false, false, false> >::value));

    have_is_nothrow_move =
      boost::is_nothrow_move_constructible<hash_nothrow_move_construct>::value;
    have_is_nothrow_move_assign =
      boost::is_nothrow_move_assignable<hash_nothrow_move_assign>::value;
    have_is_nothrow_swap =
      boost::is_nothrow_swappable<hash_nothrow_swap>::value;

// Check that the traits work when expected.
#if !defined(BOOST_NO_CXX11_NOEXCEPT) && !defined(BOOST_NO_SFINAE_EXPR) &&     \
  !BOOST_WORKAROUND(BOOST_GCC_VERSION, < 40800)
    BOOST_TEST(have_is_nothrow_move);
    BOOST_TEST(have_is_nothrow_move_assign);
#endif

#if !defined(BOOST_NO_SFINAE_EXPR) && !defined(BOOST_NO_CXX11_NOEXCEPT) &&     \
  !defined(BOOST_NO_CXX11_DECLTYPE) &&                                         \
  !defined(BOOST_NO_CXX11_FUNCTION_TEMPLATE_DEFAULT_ARGS) &&                   \
  !BOOST_WORKAROUND(BOOST_GCC_VERSION, < 40700)
    BOOST_TEST(have_is_nothrow_swap);
#endif

    BOOST_LIGHTWEIGHT_TEST_OSTREAM
      << "have_is_nothrow_move: " << have_is_nothrow_move << std::endl
      << "have_is_nothrow_swap: " << have_is_nothrow_swap << std::endl;
  }

  UNORDERED_AUTO_TEST (test_noexcept) {
    if (have_is_nothrow_move) {
#ifdef BOOST_UNORDERED_FOA_TESTS
      BOOST_TEST((boost::is_nothrow_move_constructible<
        boost::unordered_flat_set<int> >::value));
      BOOST_TEST((boost::is_nothrow_move_constructible<
        boost::unordered_flat_map<int, int> >::value));
#else
      BOOST_TEST((boost::is_nothrow_move_constructible<
        boost::unordered_set<int> >::value));
      BOOST_TEST((boost::is_nothrow_move_constructible<
        boost::unordered_multiset<int> >::value));
      BOOST_TEST((boost::is_nothrow_move_constructible<
        boost::unordered_map<int, int> >::value));
      BOOST_TEST((boost::is_nothrow_move_constructible<
        boost::unordered_multimap<int, int> >::value));
#endif
    }

#ifdef BOOST_UNORDERED_FOA_TESTS
    BOOST_TEST(
      (!boost::is_nothrow_move_constructible<
        boost::unordered_flat_set<int, hash_possible_exception> >::value));
    BOOST_TEST(
      (!boost::is_nothrow_move_constructible<boost::unordered_flat_set<int,
          boost::hash<int>, equal_to_possible_exception> >::value));
#else
    BOOST_TEST((!boost::is_nothrow_move_constructible<
                boost::unordered_set<int, hash_possible_exception> >::value));
    BOOST_TEST(
      (!boost::is_nothrow_move_constructible<boost::unordered_multiset<int,
          boost::hash<int>, equal_to_possible_exception> >::value));
#endif
  }

  UNORDERED_AUTO_TEST (test_nothrow_move_when_noexcept) {
#ifdef BOOST_UNORDERED_FOA_TESTS
    typedef boost::unordered_flat_set<int, hash_nothrow_move_construct,
      equal_to_nothrow_move_construct>
      throwing_set;
#else
    typedef boost::unordered_set<int, hash_nothrow_move_construct,
      equal_to_nothrow_move_construct>
      throwing_set;
#endif

    if (have_is_nothrow_move) {
      BOOST_TEST(boost::is_nothrow_move_constructible<throwing_set>::value);
    }

    throwing_test_exception = false;

    throwing_set x1;
    x1.insert(10);
    x1.insert(50);

    try {
      throwing_test_exception = true;

      throwing_set x2 = boost::move(x1);
      BOOST_TEST(x2.size() == 2);
      BOOST_TEST(*x2.begin() == 10 || *x2.begin() == 50);
      BOOST_TEST(have_is_nothrow_move);
    } catch (test_exception) {
      BOOST_TEST(!have_is_nothrow_move);
    }

    throwing_test_exception = false;
  }

  template <class T>
  void test_nothrow_move_assign_when_noexcept(T*, test::random_generator const&)
  {
    {
      if (have_is_nothrow_move_assign) {
        BOOST_TEST(boost::is_nothrow_move_assignable<T>::value);
      }

      throwing_test_exception = false;

      T x1;
      T x2;
      x1.insert(10);
      x1.insert(50);
      for (int i = 0; i < 100; ++i) {
        x2.insert(i);
      }

      try {
        throwing_test_exception = true;

        x2 = boost::move(x1);
        BOOST_TEST(x2.size() == 2);
        BOOST_TEST(*x2.begin() == 10 || *x2.begin() == 50);
        BOOST_TEST(have_is_nothrow_move_assign);
      } catch (test_exception) {
        BOOST_TEST(!have_is_nothrow_move_assign);
      }

      throwing_test_exception = false;
    }

    {
      if (have_is_nothrow_move_assign) {
        BOOST_TEST(boost::is_nothrow_move_assignable<T>::value);
      }

      throwing_test_exception = false;

      T x1;
      T x2;
      x1.insert(10);
      x1.insert(50);
      for (int i = 0; i < 100; ++i) {
        x2.insert(i);
      }

      try {
        throwing_test_exception = true;

        x1 = boost::move(x2);
        BOOST_TEST(x1.size() == 100);
        BOOST_TEST(have_is_nothrow_move_assign);
      } catch (test_exception) {
        BOOST_TEST(!have_is_nothrow_move_assign);
      }

      throwing_test_exception = false;
    }
  }

  UNORDERED_AUTO_TEST (test_nothrow_swap_when_noexcept) {
#if BOOST_UNORDERED_FOA_TESTS
    typedef boost::unordered_flat_set<int, hash_nothrow_swap,
      equal_to_nothrow_swap>
      throwing_set;
#else
    typedef boost::unordered_set<int, hash_nothrow_swap, equal_to_nothrow_swap>
      throwing_set;
#endif

    if (have_is_nothrow_swap) {
      BOOST_TEST(boost::is_nothrow_swappable<throwing_set>::value);
    }

    throwing_test_exception = false;

    throwing_set x1;
    throwing_set x2;
    x1.insert(10);
    x1.insert(50);
    for (int i = 0; i < 100; ++i) {
      x2.insert(i);
    }

    try {
      throwing_test_exception = true;

      x1.swap(x2);
      BOOST_TEST(x1.size() == 100);
      BOOST_TEST(x2.size() == 2);
      BOOST_TEST(*x2.begin() == 10 || *x2.begin() == 50);
      BOOST_TEST(have_is_nothrow_swap);
    } catch (test_exception) {
      BOOST_TEST(!have_is_nothrow_swap);
    }

    throwing_test_exception = false;
  }
}

#if defined(BOOST_MSVC)
#pragma warning(pop)
#endif

template <class T> class allocator1
{
  BOOST_COPYABLE_AND_MOVABLE(allocator1)
  allocator1 operator=(BOOST_COPY_ASSIGN_REF(allocator1));
  allocator1 operator=(BOOST_RV_REF(allocator1));

public:
  typedef T value_type;

  allocator1() {}
  allocator1(allocator1 const&) {}

  template <class U> allocator1(allocator1<U> const&) {}

  T* allocate(std::size_t n)
  {
    noexcept_tests::test_throw("Allocate");
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) { ::operator delete(p); }

  friend bool operator==(allocator1 const&, allocator1 const&) { return true; }
  friend bool operator!=(allocator1 const&, allocator1 const&) { return false; }
};

template <class T> class allocator2
{
  BOOST_COPYABLE_AND_MOVABLE(allocator2)
  allocator2 operator=(BOOST_COPY_ASSIGN_REF(allocator2));
public:
  typedef T value_type;
  typedef boost::true_type propagate_on_container_move_assignment;

  allocator2() {}
  allocator2(allocator2 const&) {}

  template <class U> allocator2(allocator2<U> const&) {}

  allocator2& operator=(BOOST_RV_REF(allocator2)) { return *this; }

  T* allocate(std::size_t n)
  {
    noexcept_tests::test_throw("Allocate");
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) { ::operator delete(p); }

  friend bool operator==(allocator2 const&, allocator2 const&) { return true; }
  friend bool operator!=(allocator2 const&, allocator2 const&) { return false; }
};

UNORDERED_AUTO_TEST (prelim_allocator_checks) {
  BOOST_TEST(boost::allocator_is_always_equal<allocator1<int> >::type::value);
  BOOST_TEST(!boost::allocator_propagate_on_container_move_assignment<
             allocator1<int> >::type::value);

  BOOST_TEST(boost::allocator_is_always_equal<allocator2<int> >::type::value);
  BOOST_TEST(boost::allocator_propagate_on_container_move_assignment<
             allocator2<int> >::type::value);
}

using test::default_generator;

#ifdef BOOST_UNORDERED_FOA_TESTS
boost::unordered_flat_set<int, noexcept_tests::hash_nothrow_move_assign,
  noexcept_tests::equal_to_nothrow_move_assign, allocator1<int> >*
  throwing_set_alloc1;

boost::unordered_flat_set<int, noexcept_tests::hash_nothrow_move_assign,
  noexcept_tests::equal_to_nothrow_move_assign, allocator2<int> >*
  throwing_set_alloc2;

UNORDERED_TEST(test_nothrow_move_assign_when_noexcept,
  ((throwing_set_alloc1)(throwing_set_alloc2))((default_generator)))
#else
boost::unordered_set<int, noexcept_tests::hash_nothrow_move_assign,
  noexcept_tests::equal_to_nothrow_move_assign, allocator1<int> >*
  throwing_set_alloc1;

boost::unordered_set<int, noexcept_tests::hash_nothrow_move_assign,
  noexcept_tests::equal_to_nothrow_move_assign, allocator2<int> >*
  throwing_set_alloc2;

UNORDERED_TEST(test_nothrow_move_assign_when_noexcept,
  ((throwing_set_alloc1)(throwing_set_alloc2))((default_generator)))
#endif

RUN_TESTS()
