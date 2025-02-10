// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/detached.hpp>

#include <boost/asio/steady_timer.hpp>

using namespace boost;

#include "test.hpp"

BOOST_AUTO_TEST_SUITE(detached);

cobalt::detached d(bool & done)
{
  asio::steady_timer tim{co_await cobalt::this_coro::executor, std::chrono::milliseconds(10)};
  BOOST_CHECK_NO_THROW(co_await tim.async_wait(cobalt::use_op));
  done = true;
}

BOOST_AUTO_TEST_CASE(detached)
{
  asio::io_context ctx;
  cobalt::this_thread::set_executor(ctx.get_executor());

  bool done = false;
  d(done);

  BOOST_CHECK(!done);
  ctx.run();
  BOOST_CHECK(done);
}

BOOST_AUTO_TEST_SUITE_END();