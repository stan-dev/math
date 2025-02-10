
// Copyright 2007-2009 Daniel James.
// Copyright 2022-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"
#include <string>

namespace at_tests {

  template <class X> static void at_tests(X*)
  {
    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "Create Map" << std::endl;

    X x;
    X const& x_const(x);

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

#ifdef BOOST_UNORDERED_FOA_TESTS
  static boost::unordered_flat_map<std::string, int>* test_map;
  static boost::unordered_node_map<std::string, int>* test_node_map;

  // clang-format off
  UNORDERED_TEST(at_tests, ((test_map)(test_node_map)))
  // clang-format on
#else
  static boost::unordered_map<std::string, int>* test_map;

  // clang-format off
  UNORDERED_TEST(at_tests, ((test_map)))
  // clang-format on
#endif
} // namespace at_tests

RUN_TESTS()
