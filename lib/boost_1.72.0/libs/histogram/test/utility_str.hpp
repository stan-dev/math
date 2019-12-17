// Copyright 2019 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_HISTOGRAM_TEST_UTILITY_STR_HPP
#define BOOST_HISTOGRAM_TEST_UTILITY_STR_HPP

#include <sstream>
#include <string>

template <class T>
std::string str(const T& t) {
  std::ostringstream os;
  os << t;
  return os.str();
}

#endif
