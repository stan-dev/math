// Copyright (C) 2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>

#if defined(BOOST_GCC)
// Spurious maybe-uninitialized warnings with allocators contained
// in node handles.
// Maybe related to https://gcc.gnu.org/bugzilla/show_bug.cgi?id=108230
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#ifdef BOOST_UNORDERED_CFOA_TESTS
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>
#else
#include "../helpers/unordered.hpp"
#endif

#include "../helpers/test.hpp"
#include "../helpers/replace_allocator.hpp"
  
#include <boost/config/workaround.hpp>
#include <boost/core/allocator_access.hpp>
#include <memory>
#include <type_traits>


namespace {
  template <class T> struct nonassignable_allocator
  {
    using value_type = T;

    nonassignable_allocator() = default;
    nonassignable_allocator(nonassignable_allocator const&) = default;

    template <class U>
    nonassignable_allocator(nonassignable_allocator<U> const&) {}

    nonassignable_allocator& operator=(nonassignable_allocator const&) = delete;

    T* allocate(std::size_t n)
    {
      return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t) { ::operator delete(p); }

    bool operator==(nonassignable_allocator const&) const { return true; }
    bool operator!=(nonassignable_allocator const&) const { return false; }
  };

  template <class T> struct pocx_allocator
  {
    int x_;

    using value_type = T;
    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;

    pocx_allocator() : x_{-1} {}
    pocx_allocator(pocx_allocator const&) = default;
    pocx_allocator(int const x) : x_{x} {}

    template <class U>
    pocx_allocator(pocx_allocator<U> const& rhs) : x_{rhs.x_}
    {
    }

    pocx_allocator& operator=(pocx_allocator const&) = default;

    T* allocate(std::size_t n)
    {
      return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t) { ::operator delete(p); }

    bool operator==(pocx_allocator const& rhs) const { return x_ == rhs.x_; }
    bool operator!=(pocx_allocator const& rhs) const { return x_ != rhs.x_; }
  };

  template<typename X, typename Allocator>
  void node_handle_allocator_tests(
    X*, std::pair<Allocator, Allocator> allocators)
  {
    using value_type = typename X::value_type;
    using replaced_allocator_container = test::replace_allocator<X, Allocator>;
    using node_type = typename replaced_allocator_container::node_type;

    replaced_allocator_container x1(allocators.first);
    node_type nh;

    x1.emplace(value_type());
    nh = x1.extract(0);

    BOOST_TEST(!nh.empty());
    BOOST_TEST(nh.get_allocator() == x1.get_allocator());

    replaced_allocator_container x2(allocators.second);

    x2.emplace(value_type());
    nh = x2.extract(0);

    BOOST_TEST(!nh.empty());
    BOOST_TEST(nh.get_allocator() == x2.get_allocator());
  }

  template<typename X, typename Allocator>
  void node_handle_allocator_swap_tests(
    X*, std::pair<Allocator, Allocator> allocators)
  {
    using value_type = typename X::value_type;
    using replaced_allocator_container = test::replace_allocator<X, Allocator>;
    using node_type = typename replaced_allocator_container::node_type;

    replaced_allocator_container x1(allocators.first), x2(allocators.second);
    x1.emplace(value_type());
    x2.emplace(value_type());

    node_type nh1, nh2;

    nh1 = x1.extract(0);
    swap(nh1, nh2);

    BOOST_TEST(nh1.empty());
    BOOST_TEST(!nh2.empty());
    BOOST_TEST(nh2.get_allocator() == x1.get_allocator());

    nh1 = x2.extract(0);
    swap(nh1, nh2);

    BOOST_TEST(!nh1.empty());
    BOOST_TEST(nh1.get_allocator() == x1.get_allocator());
    BOOST_TEST(!nh2.empty());
    BOOST_TEST(nh2.get_allocator() == x2.get_allocator());
  }

#if BOOST_WORKAROUND(BOOST_MSVC, <= 1900)
#pragma warning(push)
#pragma warning(disable : 4592) // symbol will be dynamically initialized
#endif

  std::pair<
    std::allocator<int>, std::allocator<int> > test_std_allocators({},{});
  std::pair<
    nonassignable_allocator<int>,
    nonassignable_allocator<int> > test_nonassignable_allocators({},{});
  std::pair<
    pocx_allocator<int>, pocx_allocator<int> > test_pocx_allocators(5,6);

#if BOOST_WORKAROUND(BOOST_MSVC, <= 1900)
#pragma warning(pop) // C4592
#endif

#if defined(BOOST_UNORDERED_FOA_TESTS)
  boost::unordered_node_map<int, int>* test_map;
  boost::unordered_node_set<int>* test_set;
#elif defined(BOOST_UNORDERED_CFOA_TESTS)
  boost::concurrent_node_map<int, int>* test_map;
  boost::concurrent_node_set<int>* test_set;
#else
  boost::unordered_map<int, int>* test_map;
  boost::unordered_set<int>* test_set;
#endif
} // namespace

// clang-format off
UNORDERED_TEST(
  node_handle_allocator_tests,
  ((test_map)(test_set))
  ((test_std_allocators)(test_nonassignable_allocators)
   (test_pocx_allocators)))

UNORDERED_TEST(
  node_handle_allocator_swap_tests,
  ((test_map)(test_set))
  ((test_std_allocators)(test_pocx_allocators)))
// clang-format on

RUN_TESTS()
