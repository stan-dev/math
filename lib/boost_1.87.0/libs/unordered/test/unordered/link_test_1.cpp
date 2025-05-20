
// Copyright 2006-2009 Daniel James.
// Copyright 2022-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#ifdef BOOST_UNORDERED_FOA_TESTS

#include <boost/unordered/concurrent_flat_map.hpp>

void foo(boost::unordered_flat_set<int>&, boost::unordered_flat_map<int, int>&,
  boost::unordered_node_set<int>&, boost::unordered_node_map<int, int>&,
  boost::concurrent_flat_map<int, int>&);

int main()
{
  boost::unordered_flat_set<int> x1;
  boost::unordered_flat_map<int, int> x2;
  boost::unordered_node_set<int> x3;
  boost::unordered_node_map<int, int> x4;
  boost::concurrent_flat_map<int, int> x5;

  foo(x1, x2, x3, x4, x5);

  return 0;
}
#else
void foo(boost::unordered_set<int>&, boost::unordered_map<int, int>&,
  boost::unordered_multiset<int>&, boost::unordered_multimap<int, int>&);

int main()
{
  boost::unordered_set<int> x1;
  boost::unordered_map<int, int> x2;
  boost::unordered_multiset<int> x3;
  boost::unordered_multimap<int, int> x4;

  foo(x1, x2, x3, x4);

  return 0;
}
#endif
