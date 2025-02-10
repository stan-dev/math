// Copyright (c) 2024 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/experimental/context.hpp>
#include <boost/cobalt/op.hpp>
#include <boost/cobalt/task.hpp>
#include <boost/cobalt/generator.hpp>

#include "../test.hpp"

using namespace boost::cobalt;

BOOST_AUTO_TEST_SUITE(context_);

BOOST_AUTO_TEST_CASE(basics)
{
  boost::cobalt::experimental::detail::context_frame<int> ff;

  auto hh = std::coroutine_handle<int>::from_address(&ff);
  BOOST_CHECK(!hh.done());
  ff.resume_ = nullptr;
  BOOST_CHECK(hh.done());

  BOOST_CHECK(&ff.promise == &hh.promise());
  BOOST_CHECK(std::coroutine_handle<int>::from_promise(ff.promise).address() == &ff);
}

int stackful_task(experimental::context<task<int>> c, int init)
{
  static_assert(experimental::detail::has_await_transform<decltype(c)::promise_type>);
  BOOST_CHECK(c.await(this_coro::executor) == c.promise().get_executor());

  bool done = false;
  boost::asio::post(c.promise().get_executor(),
                    [&]{done = true;});

  BOOST_CHECK(!done);
  c.await(boost::asio::post(use_op));
  BOOST_CHECK(done);
  return init * 2;
}

CO_TEST_CASE(task_)
{
  task<int> t = experimental::make_context(&stackful_task, 12);

  BOOST_CHECK(24 == co_await t);
}

void stackful_task_void(experimental::context<task<void>> c)
{
  static_assert(experimental::detail::has_await_transform<decltype(c)::promise_type>);
  BOOST_CHECK(c.await(this_coro::executor) == c.promise().get_executor());
}

CO_TEST_CASE(task_void)
{
  co_await experimental::make_context(&stackful_task_void);
}

int stackful_generator(experimental::context<generator<int>> c)
{
  static_assert(experimental::detail::has_await_transform<decltype(c)::promise_type>);
  BOOST_CHECK(c.await(this_coro::executor) == c.promise().get_executor());

  c.yield(1);
  c.yield(2);
  c.yield(3);
  c.yield(4);

  return 5;
}

CO_TEST_CASE(generator_)
{
  generator<int> g = experimental::make_context(&stackful_generator);
  BOOST_CHECK(1 == co_await g);
  BOOST_CHECK(2 == co_await g);
  BOOST_CHECK(3 == co_await g);
  BOOST_CHECK(4 == co_await g);
  BOOST_CHECK(5 == co_await g);
}

BOOST_AUTO_TEST_SUITE_END();
