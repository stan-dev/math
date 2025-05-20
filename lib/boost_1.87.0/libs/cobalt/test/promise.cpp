//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/promise.hpp>
#include <boost/cobalt/op.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"
#include <boost/asio/detached.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/awaitable.hpp>

#include <boost/asio/bind_cancellation_slot.hpp>

using namespace boost;

BOOST_AUTO_TEST_SUITE(promise);

cobalt::promise<void> test0()
{
    co_return;
}

cobalt::promise<double> test2(int i)
{
    co_await test0();
    co_return i;
}


cobalt::promise<int> test1(asio::any_io_executor exec)
{
    co_await test2(42);
    co_await test2(42);

    BOOST_CHECK(test2(42).get() == 42);

    co_await asio::post(exec, cobalt::use_op);
    co_return 452;
}


CO_TEST_CASE(cobalt_1)
{
    co_await test1(co_await asio::this_coro::executor);
    co_return;
}


cobalt::promise<void> should_unwind(asio::io_context & ctx)
{
  co_await asio::post(ctx, cobalt::use_op);
}

BOOST_AUTO_TEST_CASE(unwind)
{
  asio::io_context ctx;
  boost::cobalt::this_thread::set_executor(ctx.get_executor());
  +should_unwind(ctx);
}

cobalt::promise<int> return_(std::size_t ms)
{
  return cobalt::noop(1234);
}

cobalt::promise<int> return_(std::size_t ms, asio::executor_arg_t,
                            boost::cobalt::executor )
{
  co_return 1234u;
}

cobalt::promise<std::size_t> delay_r(asio::io_context &ctx, std::size_t ms,
                                    asio::executor_arg_t, asio::any_io_executor )
{
  auto tim = cobalt::use_op.as_default_on(asio::steady_timer(ctx, std::chrono::milliseconds{ms}));
  co_await tim.async_wait();
  co_return ms;
}


cobalt::promise<std::size_t> delay_r(asio::any_io_executor exec, std::size_t ms)
{
   auto tim = cobalt::use_op.as_default_on(asio::steady_timer(exec, std::chrono::milliseconds{ms}));
  co_await tim.async_wait();
  co_return ms;
}

cobalt::promise<void> throw_()
{
  throw std::runtime_error("throw_");
  co_return ;
}

cobalt::promise<void> throw_post()
{
  co_await asio::post(cobalt::this_thread::get_executor(), cobalt::use_op);
  throw std::runtime_error("throw_");
  co_return ;
}


BOOST_AUTO_TEST_CASE(get)
{
  BOOST_CHECK_THROW(throw_().get(), std::exception);
}

cobalt::promise<void> delay_v(asio::io_context &ctx, std::size_t ms)
{
  asio::steady_timer tim(ctx, std::chrono::milliseconds{ms});
  co_await tim.async_wait(boost::cobalt::use_op);
}


CO_TEST_CASE(cancel_int)
{
  BOOST_CHECK_THROW(co_await throw_(), std::exception);
}



CO_TEST_CASE(throw_cpl_delay)
{
  BOOST_CHECK_THROW(co_await throw_post(), std::exception);
}

CO_TEST_CASE(stop_)
try
{
  BOOST_CHECK_THROW(
    co_await
        []() -> cobalt::promise<void>
        {
          co_await stop();
        }(), boost::system::system_error);
}
catch(...)
{
}

struct promise_move_only
{
  promise_move_only() = default;
  promise_move_only(promise_move_only &&) = default;
  promise_move_only & operator=(promise_move_only &&) = default;
};

cobalt::promise<promise_move_only> pro_move_only_test()
{
  co_return promise_move_only{};
}

CO_TEST_CASE(move_only)
{
  auto p = pro_move_only_test();
  co_await p;
  p = pro_move_only_test();
  co_await p;
}

BOOST_AUTO_TEST_SUITE_END();