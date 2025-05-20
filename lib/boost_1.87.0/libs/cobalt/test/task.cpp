//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/task.hpp>
#include <boost/cobalt/main.hpp>
#include <boost/cobalt/op.hpp>
#include <boost/cobalt/spawn.hpp>
#include <boost/cobalt/join.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"
#include <boost/asio/cancel_after.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/bind_cancellation_slot.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/strand.hpp>
#include <boost/asio/thread_pool.hpp>

#include <boost/core/ignore_unused.hpp>

#include <new>

using namespace boost;

BOOST_AUTO_TEST_SUITE(task);

namespace
{

cobalt::task<void> test0()
{
    return cobalt::noop<void>();
}

cobalt::task<double> test2(int i)
{
    co_await cobalt::noop();
    co_await test0();
    co_return i;
}

cobalt::task<int> test1()
{
    co_await test2(42);
    co_await test2(42);


    co_await asio::post(co_await cobalt::this_coro::executor, cobalt::use_task);
    co_return 452;
}
}

BOOST_AUTO_TEST_CASE(test_1)
{

    bool done = false;

    asio::io_context ctx;
    asio::steady_timer tim{ctx};
    cobalt::this_thread::set_executor(ctx.get_executor());

    cobalt::spawn(
          ctx.get_executor(),
          test1(),
          [&](std::exception_ptr ex, int res)
          {
            BOOST_CHECK(ex == nullptr);
            BOOST_CHECK(res == 452);
            done = true;
          });

    ctx.run();
    BOOST_CHECK(done);
}

CO_TEST_CASE(cobalt_1)
{
    co_await test1();
    co_return;
}


static cobalt::task<void> should_unwind(asio::io_context & ctx)
{
  co_await asio::post(ctx, cobalt::use_op);
}

BOOST_AUTO_TEST_CASE(unwind)
{
  asio::io_context ctx;
  boost::cobalt::this_thread::set_executor(ctx.get_executor());
  boost::ignore_unused(should_unwind(ctx));
}
namespace
{


cobalt::task<int> return_([[maybe_unused]] std::size_t ms,
                         asio::executor_arg_t, boost::cobalt::executor )
{
  co_return 1234u;
}

cobalt::task<void> delay_r(asio::io_context &ctx, std::size_t ms)
{
  auto tim = cobalt::use_op.as_default_on(asio::steady_timer(ctx, std::chrono::milliseconds{ms}));
  co_await tim.async_wait();
}

cobalt::task<void> throw_()
{
  throw std::runtime_error("throw_");
  co_return ;
}

cobalt::task<void> throw_post()
{
  co_await asio::post(cobalt::this_thread::get_executor(), cobalt::use_op);
  throw std::runtime_error("throw_");
  co_return ;
}

}


BOOST_AUTO_TEST_CASE(cancel_void)
{
  asio::io_context ctx;
  cobalt::this_thread::set_executor(ctx.get_executor());
  asio::cancellation_signal signal;


  spawn(ctx, delay_r(ctx, 10000u), asio::bind_cancellation_slot(
      signal.slot(),
      [](std::exception_ptr ep)
      {
        BOOST_CHECK(ep != nullptr);
      }));

  asio::post(ctx, [&]{signal.emit(asio::cancellation_type::all);});

  spawn(ctx, return_(1234u, asio::executor_arg, ctx.get_executor()),
        [](std::exception_ptr ep, std::size_t n)
        {
          BOOST_CHECK(ep == nullptr);
          BOOST_CHECK(n == 1234u);
        });


  ctx.run();
}

static cobalt::task<void> delay_v(asio::io_context &ctx, std::size_t ms)
{
  asio::steady_timer tim(ctx, std::chrono::milliseconds{ms});
  co_await tim.async_wait(cobalt::use_op);
}


BOOST_AUTO_TEST_CASE(cancel_int)
{
  asio::io_context ctx;
  cobalt::this_thread::set_executor(ctx.get_executor());
  asio::cancellation_signal signal;

  spawn(ctx, delay_v(ctx, 10000u), asio::bind_cancellation_slot(
      signal.slot(),
      [](std::exception_ptr ep)
      {
        BOOST_CHECK(ep != nullptr);
      }));

  asio::post(ctx, [&]{signal.emit(asio::cancellation_type::all);});
  spawn(ctx, throw_(),
          [](std::exception_ptr ep)
          {
            BOOST_CHECK(ep != nullptr);
          });


  ctx.run();
}


BOOST_AUTO_TEST_CASE(throw_cpl)
{
  asio::io_context ctx;
  cobalt::this_thread::set_executor(ctx.get_executor());
  asio::cancellation_signal signal;

  spawn(ctx, throw_(),
        [](std::exception_ptr ep)
        {
          std::rethrow_exception(ep);
        });


  BOOST_CHECK_THROW(ctx.run(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(throw_cpl_delay)
{
  asio::io_context ctx;
  cobalt::this_thread::set_executor(ctx.get_executor());
  asio::cancellation_signal signal;

  spawn(ctx, throw_post(),
        [](std::exception_ptr ep)
        {
          std::rethrow_exception(ep);
        });


  BOOST_CHECK_THROW(ctx.run(), std::runtime_error);
}


CO_TEST_CASE(stop_)
{
  BOOST_CHECK_THROW(
      co_await
          []() -> cobalt::task<int>
          {
            co_await stop();
            co_return 42;
          }(), boost::system::system_error);
}


cobalt::task<void> throw_if_test(asio::cancellation_signal & sig)
{

  BOOST_CHECK(co_await cobalt::this_coro::cancelled
        == asio::cancellation_type::none);
  sig.emit(asio::cancellation_type::terminal);
  BOOST_CHECK(co_await cobalt::this_coro::cancelled
        == asio::cancellation_type::terminal);
  BOOST_CHECK_THROW(co_await asio::post(cobalt::use_op), boost::system::system_error);
}


BOOST_AUTO_TEST_CASE(throw_if_cancelled)
{
  asio::cancellation_signal sig;

  asio::io_context ctx;
  boost::cobalt::this_thread::set_executor(ctx.get_executor());
  cobalt::spawn(ctx, throw_if_test(sig),
               asio::bind_cancellation_slot(sig.slot(), asio::detached));
  ctx.run();
}

CO_TEST_CASE(reawait)
{
  auto t = test0();
  co_await std::move(t);
  BOOST_CHECK_NO_THROW(co_await std::move(t));
}


cobalt::task<int> test_strand1(asio::any_io_executor exec)
{
  BOOST_ASSERT(exec == co_await cobalt::this_coro::executor);
  co_await asio::post(co_await cobalt::this_coro::executor, cobalt::use_task);
  co_return 31;
}

cobalt::task<void> test_strand()
{
  auto e = co_await cobalt::this_coro::executor;
  co_await cobalt::join(test_strand1(e), test_strand1(e), test_strand1(e), test_strand1(e));
}

#if !defined(BOOST_COBALT_USE_IO_CONTEXT)

BOOST_AUTO_TEST_CASE(stranded)
{

  asio::thread_pool ctx;
  boost::cobalt::this_thread::set_executor(ctx.get_executor());
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  cobalt::spawn(asio::make_strand(ctx.get_executor()), test_strand(),[](std::exception_ptr ep) { if (ep) std::rethrow_exception(ep); });
  BOOST_CHECK_NO_THROW(ctx.join());
}

#endif

struct task_move_only
{
  task_move_only() = default;
  task_move_only(task_move_only &&) = default;
  task_move_only & operator=(task_move_only &&) = delete;
};

cobalt::task<task_move_only> task_move_only_test()
{
  co_return task_move_only{};
}

CO_TEST_CASE(move_only)
{
  co_await task_move_only_test();
}

#if !defined(BOOST_COBALT_USE_IO_CONTEXT) && !defined(BOOST_COBALT_CUSTOM_EXECUTOR)

BOOST_AUTO_TEST_CASE(cancel_task_)
{
  asio::thread_pool ctx{1};
  asio::cancellation_signal sig;
  cobalt::spawn(ctx, test0(), asio::bind_cancellation_slot(sig.slot(), asio::detached));
  sig.emit(asio::cancellation_type::all);
  ctx.join();
}

#endif


CO_TEST_CASE(cancel_task_after)
{
  try
  {
    auto exec = co_await cobalt::this_coro::executor;
    co_await cobalt::spawn(
        exec,
        []() -> cobalt::task<void>
        {
          asio::steady_timer t{co_await cobalt::this_coro::executor, std::chrono::steady_clock::time_point::max()};
          co_await t.async_wait();
        }(),
        asio::cancel_after(std::chrono::milliseconds(1)));

    BOOST_FAIL("Unreachable");
  }
  catch(system::system_error & se)
  {
    BOOST_CHECK(se.code() == asio::error::operation_aborted);
  }
  catch(...)
  {
    BOOST_FAIL("Unreachable");
  }
}




BOOST_AUTO_TEST_SUITE_END();