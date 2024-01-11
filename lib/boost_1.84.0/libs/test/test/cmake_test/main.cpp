//
// Copyright (c) 2019-2021 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE basic_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(some_suite)

BOOST_AUTO_TEST_CASE(case1)
{
  BOOST_TEST(true);
}

BOOST_AUTO_TEST_SUITE_END();
