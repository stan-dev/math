
// Copyright 2007-2009 Daniel James.
// Copyright 2022 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"
#include <string>

namespace at_tests {

  UNORDERED_AUTO_TEST (at_tests) {
    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Create Map" << std::endl;

#ifdef BOOST_UNORDERED_FOA_TESTS
    boost::unordered_flat_map<std::string, int> x;
    boost::unordered_flat_map<std::string, int> const& x_const(x);
#else
    boost::unordered_map<std::string, int> x;
    boost::unordered_map<std::string, int> const& x_const(x);
#endif

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Check empty container" << std::endl;

    try {
      x.at("one");
      BOOST_ERROR("Should have thrown.");
    } catch (std::out_of_range&) {
    }

    try {
      x_const.at("one");
      BOOST_ERROR("Should have thrown.");
    } catch (std::out_of_range&) {
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Add elements" << std::endl;

    x["one"] = 1;
    x["two"] = 2;

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Check existing elements" << std::endl;

    BOOST_TEST(x.at("one") == 1);
    BOOST_TEST(x.at("two") == 2);
    BOOST_TEST(x_const.at("one") == 1);
    BOOST_TEST(x_const.at("two") == 2);

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Check missing element" << std::endl;

    try {
      x.at("three");
      BOOST_ERROR("Should have thrown.");
    } catch (std::out_of_range&) {
    }

    try {
      x_const.at("three");
      BOOST_ERROR("Should have thrown.");
    } catch (std::out_of_range&) {
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Finished" << std::endl;
  }
}

RUN_TESTS()
