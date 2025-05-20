
// Copyright 2006-2009 Daniel James.
// Copyright 2022 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../objects/exception.hpp"

#ifdef BOOST_UNORDERED_FOA_TESTS
typedef boost::unordered_flat_set<test::exception::object,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator<test::exception::object> >
  test_set;

typedef boost::unordered_flat_map<test::exception::object,
  test::exception::object, test::exception::hash, test::exception::equal_to,
  test::exception::allocator2<test::exception::object> >
  test_map;

typedef boost::unordered_flat_set<
  std::pair<test::exception::object, test::exception::object>,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator<test::exception::object> >
  test_pair_set;

typedef boost::unordered_node_set<test::exception::object,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator<test::exception::object> >
  test_node_set;

typedef boost::unordered_node_map<test::exception::object,
  test::exception::object, test::exception::hash, test::exception::equal_to,
  test::exception::allocator2<test::exception::object> >
  test_node_map;

typedef boost::unordered_node_set<
  std::pair<test::exception::object, test::exception::object>,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator<test::exception::object> >
  test_pair_node_set;

#define CONTAINER_SEQ (test_set)(test_map)(test_node_set)(test_node_map)
#define CONTAINER_PAIR_SEQ (test_pair_set)(test_map)(test_pair_node_set)(test_node_map)
#else
typedef boost::unordered_set<test::exception::object, test::exception::hash,
  test::exception::equal_to,
  test::exception::allocator<test::exception::object> >
  test_set;
typedef boost::unordered_multiset<test::exception::object,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator2<test::exception::object> >
  test_multiset;
typedef boost::unordered_map<test::exception::object, test::exception::object,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator2<test::exception::object> >
  test_map;
typedef boost::unordered_multimap<test::exception::object,
  test::exception::object, test::exception::hash, test::exception::equal_to,
  test::exception::allocator<test::exception::object> >
  test_multimap;
typedef boost::unordered_set<
  std::pair<test::exception::object, test::exception::object>,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator<test::exception::object> >
  test_pair_set;
typedef boost::unordered_multiset<
  std::pair<test::exception::object, test::exception::object>,
  test::exception::hash, test::exception::equal_to,
  test::exception::allocator2<test::exception::object> >
  test_pair_multiset;

#define CONTAINER_SEQ (test_set)(test_multiset)(test_map)(test_multimap)
#define CONTAINER_PAIR_SEQ                                                     \
  (test_pair_set)(test_pair_multiset)(test_map)(test_multimap)

#endif
