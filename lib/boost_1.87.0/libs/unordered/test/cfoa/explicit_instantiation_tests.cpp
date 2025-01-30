// Copyright 2024 Braden Ganetsky
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

template class boost::concurrent_flat_map<int, int>;
template class boost::concurrent_flat_set<int>;
template class boost::concurrent_node_map<int, int>;
template class boost::concurrent_node_set<int>;

int main() { return 0; }
