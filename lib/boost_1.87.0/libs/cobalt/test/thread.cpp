//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/thread.hpp>
#include <boost/asio/steady_timer.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"
#include "boost/cobalt/spawn.hpp"

boost::cobalt::thread thr()
{
  boost::asio::steady_timer tim{co_await boost::asio::this_coro::executor, std::chrono::milliseconds(100)};

  auto exec = co_await boost::asio::this_coro::executor;
  co_await tim.async_wait(boost::cobalt::use_op);
}

BOOST_AUTO_TEST_SUITE(thread);

BOOST_AUTO_TEST_CASE(run)
{
  auto t = thr();

  t.join();
  BOOST_CHECK_THROW(t.get_executor(), boost::asio::execution::bad_executor);
}


boost::cobalt::thread thr_stop()
{
  boost::asio::steady_timer tim{co_await boost::asio::this_coro::executor, std::chrono::milliseconds(100)};

#if !defined(BOOST_COBALT_USE_IO_CONTEXT)
  auto exec = co_await boost::asio::this_coro::executor;
  auto execp = exec.target<boost::asio::io_context::executor_type>();
  BOOST_ASSERT(execp != nullptr);
  auto exc = *execp;
#else
  auto exc = co_await boost::asio::this_coro::executor;
#endif

  exc.context().stop();
  co_await tim.async_wait(boost::cobalt::use_op);
}

BOOST_AUTO_TEST_CASE(stop)
{
  auto t = thr_stop();
  t.join();
}

CO_TEST_CASE(await_thread)
{
  co_await thr();

  auto th = thr_stop();
  boost::asio::steady_timer tim{co_await boost::asio::this_coro::executor, std::chrono::milliseconds(200)};
  co_await tim.async_wait(boost::cobalt::use_op);
  co_await th;
  try {
    BOOST_CHECK_THROW(co_await th, boost::system::system_error);
  } catch(...) {}
}

boost::cobalt::task<std::thread::id> on_thread()
{
  co_return std::this_thread::get_id();
}

CO_TEST_CASE(spawn_onto_thread)
{
  using namespace boost;
  auto t = thr();

  auto id = co_await cobalt::spawn(t.get_executor(), on_thread(), cobalt::use_op);
  auto id2 = t.get_id();
  auto id3 = std::this_thread::get_id();

  BOOST_CHECK(id == id2);
  BOOST_CHECK(id3 != id);

  if (t.joinable())
    t.join();
}


BOOST_AUTO_TEST_SUITE_END();