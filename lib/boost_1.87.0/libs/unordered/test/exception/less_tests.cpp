
// Copyright 2023 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "./containers.hpp"

#include "../helpers/helpers.hpp"
#include "../helpers/invariants.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/strong.hpp"
#include "../helpers/tracker.hpp"

#include <vector>

UNORDERED_AUTO_TEST (less_osx_regression) {
  DISABLE_EXCEPTIONS;
  typedef test_pair_set::value_type value_type;
  typedef test::exception::object object;

  std::vector<value_type> v;
  v.push_back(value_type(object(12, 98), object(88, 13)));
  v.push_back(value_type(object(24, 71), object(62, 84)));
  v.push_back(value_type(object(30, 0), object(5, 73)));
  v.push_back(value_type(object(34, 64), object(79, 58)));
  v.push_back(value_type(object(36, 95), object(64, 23)));
  v.push_back(value_type(object(42, 89), object(68, 44)));
  v.push_back(value_type(object(42, 26), object(93, 64)));
  v.push_back(value_type(object(86, 86), object(16, 62)));
  v.push_back(value_type(object(86, 86), object(75, 23)));
  v.push_back(value_type(object(92, 37), object(41, 90)));

  BOOST_TEST_EQ(v.size(), 10u);

  std::set<value_type, test::exception::less> s;
  s.insert(v.begin(), v.end());
  BOOST_TEST_EQ(s.size(), v.size());

  test::ordered<test_pair_set> tracker;
  test_pair_set x;
  for (std::vector<value_type>::iterator it = v.begin(); it != v.end();
       ++it) {
    x.insert(*it);
  }

  tracker.insert(v.begin(), v.end());
  tracker.compare(x);
}

RUN_TESTS()
