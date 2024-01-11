//
// Copyright (c) 2023 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/detail/monotonic_resource.hpp>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(monotonic_resource);

BOOST_AUTO_TEST_CASE(basic)
{
  char buf[1024];
  boost::cobalt::detail::monotonic_resource res{buf, sizeof(buf)};

  {
    std::vector<int, boost::cobalt::detail::monotonic_allocator<int>> vec{};

    for (int i = 0u; i < 4000; i++)
      vec.push_back(i);
  }
  {
    std::vector<int, boost::cobalt::detail::monotonic_allocator<int>> vec{&res};

    for (int i = 0u; i < 4000; i++)
      vec.push_back(i);
  }
}

BOOST_AUTO_TEST_CASE(too_small)
{
  char buf[1];
  boost::cobalt::detail::monotonic_resource res{buf, sizeof(buf)};

  {
    std::vector<int, boost::cobalt::detail::monotonic_allocator<int>> vec{};

    for (int i = 0u; i < 4000; i++)
      vec.push_back(i);
  }
  {
    std::vector<int, boost::cobalt::detail::monotonic_allocator<int>> vec{&res};

    for (int i = 0u; i < 4000; i++)
      vec.push_back(i);
  }
}

BOOST_AUTO_TEST_CASE(no_buf)
{
  boost::cobalt::detail::monotonic_resource res;

  {
    std::vector<int, boost::cobalt::detail::monotonic_allocator<int>> vec{};

    for (int i = 0u; i < 4000; i++)
      vec.push_back(i);
  }
  {
    std::vector<int, boost::cobalt::detail::monotonic_allocator<int>> vec{&res};

    for (int i = 0u; i < 4000; i++)
      vec.push_back(i);
  }
}

BOOST_AUTO_TEST_SUITE_END();