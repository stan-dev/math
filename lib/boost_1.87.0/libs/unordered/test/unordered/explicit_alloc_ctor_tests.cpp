// Copyright 2024 Joaquin M Lopez Munoz.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_UNORDERED_CFOA_TESTS
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>
#else
#include "../helpers/unordered.hpp"
#endif

#include "../helpers/helpers.hpp"
#include "../helpers/test.hpp"
#include <boost/type_traits/make_void.hpp>
#include <memory>
#include <vector>

template <class T> struct explicit_allocator
{
  typedef T value_type;

  explicit_allocator() {}

  template <class U> 
  explicit explicit_allocator(const explicit_allocator<U>&) {}

  template<class U>
  bool operator==(explicit_allocator<U> const &) const noexcept
  {
    return true;
  }

  template<class U>
  bool operator!=(explicit_allocator<U> const &) const noexcept
  {
    return false;
  }
    
  T* allocate(std::size_t n)
  {
    return std::allocator<T>().allocate(n); 
  }

  void deallocate(T* p, std::size_t n)
  {
    std::allocator<T>().deallocate(p, n); 
  }
};

template<typename Container,typename = void> struct has_extract:
std::false_type {};

template<typename Container>
struct has_extract<
  Container, boost::void_t<typename Container::insert_return_type>
>: std::true_type {};

template <class Container> void test_explicit_alloc_ctor_extract(
  Container&, std::false_type)
{
}

template <class Container> void test_explicit_alloc_ctor_extract(
  Container& c, std::true_type)
{
#ifdef BOOST_UNORDERED_CFOA_TESTS
  typename Container::key_type k;
  c.cvisit_while([&](typename Container::value_type const & x) {
    k = test::get_key<Container>(x);
    return false;
  });
  auto n = c.extract(k);
#else
  auto n = c.extract(c.begin());
#endif
  c.insert(std::move(n));
  n = c.extract(typename Container::key_type());
  c.insert(std::move(n));
}

template <class Container> void test_explicit_alloc_ctor_extract(Container& c)
{
  test_explicit_alloc_ctor_extract(c, has_extract<Container>());
}

template <class Container> void test_explicit_alloc_ctor()
{
  using key_type = typename Container::key_type;
  using value_type = typename Container::value_type;
  using hasher = typename Container::hasher;
  using allocator_type = typename Container::allocator_type;

  allocator_type                    al;
  std::initializer_list<value_type> il{};
  hasher                            h;
  std::vector<value_type>           v;
  Container                         c,
                                    c2(0),
                                    c3(v.begin(), v.end()),
                                    c4(c3),
                                    c5(std::move(c4)),
                                    c6(v.begin(), v.end(), al),
                                    c7(al),
                                    c8(c7, al),
                                    c9(std::move(c8), al),
                                    c10(il),
                                    c11(0, al),
                                    c12(0, h, al),
                                    c13(v.begin(), v.end(), 0, al),
                                    c14(v.begin(), v.end(), 0, h, al),
                                    c15(il, al),
                                    c16(il, 0, al),
                                    c17(il, 0, h, al);
  value_type                        x{};
  key_type                          k{};

  c = c2;
  c = std::move(c2);
  c = il;
  c.swap(c3);
  c.insert(x);
#ifdef BOOST_UNORDERED_CFOA_TESTS
  (void) c.visit(k, [](const value_type&) {});
#else
  (void) c.find(k);
#endif
  test_explicit_alloc_ctor_extract(c);
}

UNORDERED_AUTO_TEST (explicit_alloc_ctor) {
#if defined(BOOST_UNORDERED_CFOA_TESTS)
  test_explicit_alloc_ctor<boost::concurrent_flat_map<int, int,
    boost::hash<int>, std::equal_to<int>,
    explicit_allocator<std::pair<const int, int> > > >();
  test_explicit_alloc_ctor<boost::concurrent_node_map<int, int,
    boost::hash<int>, std::equal_to<int>,
    explicit_allocator<std::pair<const int, int> > > >();
  test_explicit_alloc_ctor<boost::concurrent_flat_set<
    int, boost::hash<int>, std::equal_to<int>, explicit_allocator<int> > >();
  test_explicit_alloc_ctor<boost::concurrent_node_set<
    int, boost::hash<int>, std::equal_to<int>, explicit_allocator<int> > >();
#elif defined(BOOST_UNORDERED_FOA_TESTS)
  test_explicit_alloc_ctor<boost::unordered_flat_map<int, int,
    boost::hash<int>, std::equal_to<int>,
    explicit_allocator<std::pair<const int, int> > > >();
  test_explicit_alloc_ctor<boost::unordered_flat_set<
    int, boost::hash<int>, std::equal_to<int>, explicit_allocator<int> > >();
  test_explicit_alloc_ctor<boost::unordered_node_map<int, int,
    boost::hash<int>, std::equal_to<int>,
    explicit_allocator<std::pair<const int, int> > > >();
  test_explicit_alloc_ctor<boost::unordered_node_set<
    int, boost::hash<int>, std::equal_to<int>, explicit_allocator<int> > >();
#else
  test_explicit_alloc_ctor<boost::unordered_map<int, int,
    boost::hash<int>, std::equal_to<int>,
    explicit_allocator<std::pair<const int, int> > > >();
  test_explicit_alloc_ctor<boost::unordered_multimap<int, int,
    boost::hash<int>, std::equal_to<int>,
    explicit_allocator<std::pair<const int, int> > > >();
  test_explicit_alloc_ctor<boost::unordered_set<
    int, boost::hash<int>, std::equal_to<int>, explicit_allocator<int> > >();
  test_explicit_alloc_ctor<boost::unordered_multiset<
    int, boost::hash<int>, std::equal_to<int>, explicit_allocator<int> > >();
#endif
}

RUN_TESTS()
