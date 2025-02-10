// Copyright 2024 Braden Ganetsky
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifdef BOOST_UNORDERED_FOA_TESTS

#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_node_map.hpp>
#include <boost/unordered/unordered_node_set.hpp>

template class boost::unordered_flat_map<int, int>;
template class boost::unordered_flat_set<int>;
template class boost::unordered_node_map<int, int>;
template class boost::unordered_node_set<int>;

#else

#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>

template class boost::unordered_map<int, int>;
template class boost::unordered_multimap<int, int>;
template class boost::unordered_multiset<int>;
template class boost::unordered_set<int>;

#endif // BOOST_UNORDERED_FOA_TESTS

int main() { return 0; }
