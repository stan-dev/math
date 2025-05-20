// Copyright 2021-2022 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"

#include <boost/config.hpp>
#include <string>

#define UNORDERED_LVALUE_QUAL &

namespace test {
  struct is_even
  {
    bool operator()(std::pair<std::string const, int>& key_value)
    {
      int const v = key_value.second;
      return (v % 2 == 0);
    }

    bool operator()(int const& value)
    {
      int const v = value;
      return (v % 2 == 0);
    }
  };

  struct is_too_large
  {
    bool operator()(std::pair<std::string const, int>& key_value)
    {
      int const v = key_value.second;
      return v >= 1000;
    }

    bool operator()(int const& value)
    {
      int const v = value;
      return v >= 1000;
    }
  };

} // namespace test

template <class UnorderedMap> void test_map_erase_if()
{
  typedef UnorderedMap map_type;
  typedef typename map_type::size_type size_type;

  map_type map;
  size_type num_erased = erase_if(map, test::is_even());
  BOOST_TEST(map.empty());
  BOOST_TEST_EQ(num_erased, 0u);

  map.emplace("a", 1);
  map.emplace("b", 2);
  map.emplace("b", 4);
  map.emplace("b", 8);
  map.emplace("b", 16);
  map.emplace("c", 3);

  size_type size = map.size();

  num_erased = erase_if(map, test::is_too_large());

  BOOST_TEST_EQ(map.size(), size);
  BOOST_TEST_EQ(num_erased, 0u);

  num_erased = erase_if(map, test::is_even());
  BOOST_TEST_EQ(map.size(), 2u);
  BOOST_TEST_EQ(num_erased, size - map.size());
}

template <class UnorderedSet> void test_set_erase_if()
{
  typedef UnorderedSet set_type;
  typedef typename set_type::size_type size_type;

  set_type set;
  size_type num_erased = erase_if(set, test::is_even());
  BOOST_TEST(set.empty());
  BOOST_TEST_EQ(num_erased, 0u);

  set.emplace(1);
  set.emplace(2);
  set.emplace(2);
  set.emplace(2);
  set.emplace(2);
  set.emplace(3);

  size_type size = set.size();

  num_erased = erase_if(set, test::is_too_large());

  BOOST_TEST_EQ(set.size(), size);
  BOOST_TEST_EQ(num_erased, 0u);

  num_erased = erase_if(set, test::is_even());
  BOOST_TEST_EQ(set.size(), 2u);
  BOOST_TEST_EQ(num_erased, size - set.size());
}

UNORDERED_AUTO_TEST (unordered_erase_if) {
#ifdef BOOST_UNORDERED_FOA_TESTS
  test_map_erase_if<boost::unordered_flat_map<std::string, int> >();
  test_set_erase_if<boost::unordered_flat_set<int> >();
  test_map_erase_if<boost::unordered_node_map<std::string, int> >();
  test_set_erase_if<boost::unordered_node_set<int> >();
#else
  test_map_erase_if<boost::unordered_map<std::string, int> >();
  test_map_erase_if<boost::unordered_multimap<std::string, int> >();

  test_set_erase_if<boost::unordered_set<int> >();
  test_set_erase_if<boost::unordered_multiset<int> >();
#endif
}

RUN_TESTS()
