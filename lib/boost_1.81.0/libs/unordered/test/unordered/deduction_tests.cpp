
// Copyright 2017-2018 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test_trait.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <iostream>
#include <vector>

#if BOOST_UNORDERED_TEMPLATE_DEDUCTION_GUIDES

#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>

struct hash_equals
{
  template <typename T> bool operator()(T const& x) const
  {
    boost::hash<T> hf;
    return hf(x);
  }

  template <typename T> bool operator()(T const& x, T const& y) const
  {
    std::equal_to<T> eq;
    return eq(x, y);
  }
};

template <typename T> struct test_allocator
{
  typedef T value_type;
  test_allocator() = default;
  template <typename T2> test_allocator(test_allocator<T2> const&) {}
  T* allocate(std::size_t n) const { return (T*)malloc(sizeof(T) * n); }
  void deallocate(T* ptr, std::size_t) const { free(ptr); }
  bool operator==(test_allocator const&) const { return true; }
  bool operator!=(test_allocator const&) const { return false; }
};

template <template <class...> class UnorderedMap> void map_tests()
{
  std::vector<std::pair<int, int> > x;
  x.push_back(std::make_pair(1, 3));
  x.push_back(std::make_pair(5, 10));
  test_allocator<std::pair<const int, int> > pair_allocator;
  hash_equals f;

  /*
   template<class InputIterator,
           class Hash = hash<iter_key_t<InputIterator>>,
           class Pred = equal_to<iter_key_t<InputIterator>>,
           class Allocator = allocator<iter_to_alloc_t<InputIterator>>>
    unordered_map(InputIterator, InputIterator, typename see below::size_type =
   see below,
                  Hash = Hash(), Pred = Pred(), Allocator = Allocator())
      -> unordered_map<iter_key_t<InputIterator>, iter_val_t<InputIterator>,
   Hash, Pred,
                       Allocator>;
  */

  {
    UnorderedMap m(x.begin(), x.end());
    BOOST_TEST_TRAIT_SAME(decltype(m), UnorderedMap<int, int>);
  }

  {
    UnorderedMap m(x.begin(), x.end(), 0, std::hash<int>());
    BOOST_TEST_TRAIT_SAME(decltype(m), UnorderedMap<int, int, std::hash<int> >);
  }

  {
    UnorderedMap m(
      x.begin(), x.end(), 0, std::hash<int>(), std::equal_to<int>());
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, std::hash<int>, std::equal_to<int> >);
  }

  {
    UnorderedMap m(x.begin(), x.end(), 0, std::hash<int>(),
      std::equal_to<int>(), pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, std::hash<int>, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  /*
  template<class Key, class T, class Hash = hash<Key>,
           class Pred = equal_to<Key>, class Allocator = allocator<pair<const
  Key, T>>>
    unordered_map(initializer_list<pair<const Key, T>>,
                  typename see below::size_type = see below, Hash = Hash(),
                  Pred = Pred(), Allocator = Allocator())
      -> unordered_map<Key, T, Hash, Pred, Allocator>;
  */

  {
    UnorderedMap m({std::pair<int const, int>(1, 2)});
    BOOST_TEST_TRAIT_SAME(decltype(m), UnorderedMap<int, int>);
  }

  {
    UnorderedMap m({std::pair<int, int>(1, 2)});
    BOOST_TEST_TRAIT_SAME(decltype(m), UnorderedMap<int, int>);
  }

  {
    UnorderedMap m({std::pair<int const, int>(1, 2)}, 0, std::hash<int>());
    BOOST_TEST_TRAIT_SAME(decltype(m), UnorderedMap<int, int, std::hash<int> >);
  }

  {
    UnorderedMap m({std::pair<int, int>(1, 2)}, 0, std::hash<int>());
    BOOST_TEST_TRAIT_SAME(decltype(m), UnorderedMap<int, int, std::hash<int> >);
  }

  {
    UnorderedMap m({std::pair<int const, int>(1, 2)}, 0, std::hash<int>(),
      std::equal_to<int>());
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, std::hash<int>, std::equal_to<int> >);
  }

  {
    UnorderedMap m(
      {std::pair<int, int>(1, 2)}, 0, std::hash<int>(), std::equal_to<int>());
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, std::hash<int>, std::equal_to<int> >);
  }

  {
    UnorderedMap m({std::pair<int const, int>(1, 2)}, 0, f, f, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, hash_equals, hash_equals,
                     test_allocator<std::pair<const int, int> > >);
  }

  {
    UnorderedMap m({std::pair<int, int>(1, 2)}, 0, f, f, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, hash_equals, hash_equals,
                     test_allocator<std::pair<const int, int> > >);
  }

  /*
  template<class InputIterator, class Allocator>
    unordered_map(InputIterator, InputIterator, typename see below::size_type,
  Allocator)
      -> unordered_map<iter_key_t<InputIterator>, iter_val_t<InputIterator>,
                       hash<iter_key_t<InputIterator>>,
  equal_to<iter_key_t<InputIterator>>,
                       Allocator>;
  */

  {
    UnorderedMap m(x.begin(), x.end(), 0u, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  /*
    template<class InputIterator, class Allocator>
    unordered_map(InputIterator, InputIterator, Allocator)
      -> unordered_map<iter_key_t<InputIterator>, iter_val_t<InputIterator>,
                       hash<iter_key_t<InputIterator>>,
    equal_to<iter_key_t<InputIterator>>,
                       Allocator>;
  */

  {
    UnorderedMap m(x.begin(), x.end(), pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  /*
  template<class InputIterator, class Hash, class Allocator>
    unordered_map(InputIterator, InputIterator, typename see below::size_type,
  Hash, Allocator)
      -> unordered_map<iter_key_t<InputIterator>, iter_val_t<InputIterator>,
  Hash,
                       equal_to<iter_key_t<InputIterator>>, Allocator>;
  */

  {
    UnorderedMap m(x.begin(), x.end(), 0u, f, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, hash_equals, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  /*
    template<class Key, class T, typename Allocator>
    unordered_map(initializer_list<pair<const Key, T>>, typename see
    below::size_type,
                  Allocator)
      -> unordered_map<Key, T, hash<Key>, equal_to<Key>, Allocator>;
  */

  {
    UnorderedMap m({std::pair<int const, int>(1, 2)}, 0, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  {
    UnorderedMap m({std::pair<int, int>(1, 2)}, 0, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  /*
  template<class Key, class T, typename Allocator>
    unordered_map(initializer_list<pair<const Key, T>>, Allocator)
      -> unordered_map<Key, T, hash<Key>, equal_to<Key>, Allocator>;
  */

  {
    UnorderedMap m({std::pair<int const, int>(1, 2)}, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  {
    UnorderedMap m({std::pair<int, int>(1, 2)}, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  /*
  template<class Key, class T, class Hash, class Allocator>
    unordered_map(initializer_list<pair<const Key, T>>, typename see
  below::size_type, Hash,
                  Allocator)
      -> unordered_map<Key, T, Hash, equal_to<Key>, Allocator>;
  */

  {
    UnorderedMap m({std::pair<int const, int>(1, 2)}, 0, f, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, hash_equals, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }

  {
    UnorderedMap m({std::pair<int, int>(1, 2)}, 0, f, pair_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(m), UnorderedMap<int, int, hash_equals, std::equal_to<int>,
                     test_allocator<std::pair<const int, int> > >);
  }
}

template <template <class...> class UnorderedSet> void set_tests()
{
  std::vector<int> y;
  y.push_back(1);
  y.push_back(2);

  hash_equals f;

  test_allocator<int> int_allocator;

  /*   template<class InputIt,
           class Hash = std::hash<typename
  std::iterator_traits<InputIt>::value_type>, class Pred =
    std::equal_to<typename std::iterator_traits<InputIt>::value_type>, class
    Alloc = std::allocator<typename std::iterator_traits<InputIt>::value_type>>
  unordered_set(InputIt, InputIt,
            typename see below ::size_type = see below,
           Hash = Hash(), Pred = Pred(), Alloc = Alloc())
    -> unordered_set<typename std::iterator_traits<InputIt>::value_type, Hash,
    Pred, Alloc>; */

  {
    UnorderedSet s(y.begin(), y.end());
    BOOST_TEST_TRAIT_SAME(decltype(s), UnorderedSet<int>);
  }

  {
    UnorderedSet s(y.begin(), y.end(), 0, std::hash<int>());
    BOOST_TEST_TRAIT_SAME(decltype(s), UnorderedSet<int, std::hash<int> >);
  }

  {
    UnorderedSet s(
      y.begin(), y.end(), 0, std::hash<int>(), std::equal_to<int>());
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, std::hash<int>, std::equal_to<int> >);
  }

  {
    UnorderedSet s(y.begin(), y.end(), 0, std::hash<int>(),
      std::equal_to<int>(), int_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, std::hash<int>, std::equal_to<int>,
                     test_allocator<int> >);
  }

  /* template<class T,
           class Hash = std::hash<T>,
           class Pred = std::equal_to<T>,
           class Alloc = std::allocator<T>>
  unordered_set(std::initializer_list<T>,
           typename see below::size_type = see below,
           Hash = Hash(), Pred = Pred(), Alloc = Alloc())
    -> unordered_set<T, Hash, Pred, Alloc>; */

  {
    UnorderedSet s({1, 2});
    BOOST_TEST_TRAIT_SAME(decltype(s), UnorderedSet<int>);
  }

  {
    UnorderedSet s({1, 2}, 0, std::hash<int>());
    BOOST_TEST_TRAIT_SAME(decltype(s), UnorderedSet<int, std::hash<int> >);
  }

  {
    UnorderedSet s({1, 2}, 0, std::hash<int>(), std::equal_to<int>());
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, std::hash<int>, std::equal_to<int> >);
  }

  {
    UnorderedSet s(
      {1, 2}, 0, std::hash<int>(), std::equal_to<int>(), int_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, std::hash<int>, std::equal_to<int>,
                     test_allocator<int> >);
  }

  {
    UnorderedSet s({1, 2}, 0, f, f, int_allocator);
    BOOST_TEST_TRAIT_SAME(decltype(s),
      UnorderedSet<int, hash_equals, hash_equals, test_allocator<int> >);
  }

  /* template<class InputIt, class Alloc>
  unordered_set(InputIt, InputIt, typename see below::size_type, Alloc)
    -> unordered_set<typename std::iterator_traits<InputIt>::value_type,
                std::hash<typename std::iterator_traits<InputIt>::value_type>,
                std::equal_to<typename
  std::iterator_traits<InputIt>::value_type>, Alloc>; */

  {
    UnorderedSet s(y.begin(), y.end(), 0u, int_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<int> >);
  }

  /* template<class InputIt, class Hash, class Alloc>
  unordered_set(InputIt, InputIt, typename see below::size_type, Hash, Alloc)
    -> unordered_set<typename std::iterator_traits<InputIt>::value_type, Hash,
               std::equal_to<typename
  std::iterator_traits<InputIt>::value_type>, Allocator>; */

  {
    UnorderedSet s(y.begin(), y.end(), 0u, f, int_allocator);
    BOOST_TEST_TRAIT_SAME(decltype(s),
      UnorderedSet<int, hash_equals, std::equal_to<int>, test_allocator<int> >);
  }

  /*   template<class T, class Allocator>
  unordered_set(std::initializer_list<T>, typename see below::size_type,
  Allocator)
    -> unordered_set<T, std::hash<T>, std::equal_to<T>, Alloc>; */

  {
    UnorderedSet s({1, 2}, 0u, int_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<int> >);
  }

  /* template<class T, class Hash, class Alloc>
  unordered_set(std::initializer_list<T>, typename see below::size_type, Hash,
  Alloc)
    -> unordered_set<T, Hash, std::equal_to<T>, Alloc>; */
  {
    UnorderedSet s({1, 2}, 0u, f, int_allocator);
    BOOST_TEST_TRAIT_SAME(decltype(s),
      UnorderedSet<int, hash_equals, std::equal_to<int>, test_allocator<int> >);
  }

  /*
    template<class InputIterator, class Allocator>
      unordered_set(InputIterator, InputIterator, Allocator)
        -> unordered_set<iter-value-type<InputIterator>,
                         hash<iter-value-type<InputIterator>>,
                         equal_to<iter-value-type<InputIterator>>,
                         Allocator>;
   */

  {
    UnorderedSet s(y.begin(), y.end(), int_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<int> >);
  }

  /*
  template<class T, class Allocator>
    unordered_set(initializer_list<T>, Allocator)
      -> unordered_set<T, hash<T>, equal_to<T>, Allocator>;
 */

  {
    UnorderedSet s({1, 2}, int_allocator);
    BOOST_TEST_TRAIT_SAME(
      decltype(s), UnorderedSet<int, boost::hash<int>, std::equal_to<int>,
                     test_allocator<int> >);
  }
}

#endif

int main()
{
  std::cout << "BOOST_UNORDERED_TEMPLATE_DEDUCTION_GUIDES: "
            << BOOST_UNORDERED_TEMPLATE_DEDUCTION_GUIDES << std::endl;

#if BOOST_UNORDERED_TEMPLATE_DEDUCTION_GUIDES
  map_tests<boost::unordered_map>();
  map_tests<boost::unordered_multimap>();
  map_tests<boost::unordered_flat_map>();
  set_tests<boost::unordered_set>();
  set_tests<boost::unordered_multiset>();
  set_tests<boost::unordered_flat_set>();

  return boost::report_errors();
#endif
}
