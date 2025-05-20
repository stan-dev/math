//
// Copyright (c) 2023 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/asio/spawn.hpp>
#include <boost/cobalt/experimental/yield_context.hpp>
#include <boost/cobalt/op.hpp>
#include <boost/cobalt/task.hpp>

#include "../test.hpp"

using namespace boost::cobalt;

BOOST_AUTO_TEST_SUITE(asio_yield_context);

struct awaitable
{
  bool await_ready() { return ready;}

  bool await_suspend(std::coroutine_handle<void> h) { this->h = h; return suspend;}

  int await_resume() {return value;}

  std::coroutine_handle<void> h;
  bool ready{false};
  bool suspend{true};
  int value;
};

awaitable aw;

int test_impl(boost::asio::yield_context ctx)
{
  return experimental::await(aw, ctx);
}

BOOST_AUTO_TEST_CASE(ready)
{
  boost::asio::io_context ctx;

  boost::asio::spawn(ctx, &test_impl,
                     [](std::exception_ptr ep, int i)
                     {
                        BOOST_CHECK_EQUAL(i, 42);
                        BOOST_CHECK(!ep);
                     });
  aw.ready = true;
  aw.suspend = false;
  aw.value = 42;

  ctx.run();

  aw.h = {};
}

BOOST_AUTO_TEST_CASE(no_suspend)
{
  boost::asio::io_context ctx;

  boost::asio::spawn(ctx, &test_impl,
                     [](std::exception_ptr ep, int i)
                     {
                       BOOST_CHECK_EQUAL(i, 43);
                       BOOST_CHECK(!ep);
                     });
  aw.ready = false;
  aw.suspend = false;
  aw.value = 43;

  ctx.run();

  aw.h = {};
}


task<void> t()
{
  co_return;
}

struct dummy_aw
{
  bool await_ready() { return false;}

  std::coroutine_handle<void> await_suspend(std::coroutine_handle<void> h) { return h;}
  void await_resume() {}
};

BOOST_AUTO_TEST_CASE(await)
{
  boost::asio::io_context ioc;

  boost::asio::spawn(ioc,
                     [&](boost::asio::basic_yield_context<boost::cobalt::executor> ctx)
                     {
                       experimental::await(dummy_aw(), ctx);
                       experimental::await(t(), ctx);
                     },
                     [](std::exception_ptr ep)
                     {
                       BOOST_CHECK(!ep);
                     });

  ioc.run();

}


BOOST_AUTO_TEST_CASE(destroy)
{
  boost::asio::io_context ctx;

  boost::asio::spawn(ctx, &test_impl,
                     [](std::exception_ptr ep, int i)
                     {
                       BOOST_CHECK(false);
                     });
  aw.ready = false;
  aw.suspend = true;
  aw.h = nullptr;
  ctx.run();
  // ASAN will complain if the yield_context doesn't get freed.
  BOOST_CHECK(aw.h != nullptr);
  aw.h.destroy();
}


BOOST_AUTO_TEST_SUITE_END();
