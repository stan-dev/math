//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/wait_group.hpp>
#include <boost/cobalt/promise.hpp>

#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/steady_timer.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"

using namespace boost;

cobalt::promise<void> gdelay(asio::any_io_executor exec,
                            std::chrono::milliseconds ms = std::chrono::milliseconds(25))
{
  if (ms.count() == 0u)
    co_return;
  if (ms == std::chrono::milliseconds ::max())
    throw std::runtime_error("wdummy_throw");

  asio::steady_timer tim{exec, ms};
  co_await tim.async_wait(cobalt::use_op);
}


BOOST_AUTO_TEST_SUITE(wait_group);

CO_TEST_CASE(grp)
{
  auto e = co_await cobalt::this_coro::executor;

  using namespace std;

  cobalt::wait_group wg;
  wg.push_back(gdelay(e));
  wg.push_back(gdelay(e, 0ms));
  wg.push_back(gdelay(e, 10ms));
  wg.push_back(gdelay(e, 20ms));

  co_await asio::post(e, cobalt::use_op);
  BOOST_CHECK(wg.size() == 4u);
  BOOST_CHECK(wg.reap() == 1u);
  BOOST_CHECK(wg.size() == 3u);

  co_await wg.wait_one();
  BOOST_CHECK(wg.size() == 2u);

  wg.cancel();
  co_await wg;

  BOOST_CHECK(wg.size() == 0u);
}


BOOST_AUTO_TEST_SUITE_END();

