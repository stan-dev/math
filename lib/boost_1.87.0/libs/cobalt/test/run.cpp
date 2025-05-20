// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/run.hpp>

#include <boost/test/unit_test.hpp>

using namespace boost;

BOOST_AUTO_TEST_SUITE(run);

BOOST_AUTO_TEST_CASE(run)
{
  BOOST_CHECK(42 == cobalt::run([]() -> cobalt::task<int> {co_return 42;}()));

  asio::io_context ctx;
  cobalt::this_thread::set_executor(ctx.get_executor());
  BOOST_CHECK(42 == cobalt::run([]() -> cobalt::task<int> {co_return 42;}()));

  BOOST_CHECK(ctx.get_executor() == cobalt::this_thread::get_executor());
}

BOOST_AUTO_TEST_SUITE_END();
