//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/op.hpp>
#include <boost/cobalt/spawn.hpp>
#include <boost/cobalt/promise.hpp>

#include <boost/asio/as_tuple.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/experimental/channel.hpp>
#include <boost/asio/steady_timer.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"

using namespace boost;

BOOST_AUTO_TEST_SUITE(op);

template<typename Timer>
struct test_wait_op : cobalt::op<system::error_code>
{
  Timer & tim;

  test_wait_op(Timer & tim) : tim(tim) {}

  void ready(cobalt::handler<system::error_code> h)
  {
    if (tim.expiry() < Timer::clock_type::now())
      h({});
  }
  void initiate(cobalt::completion_handler<system::error_code> complete)
  {
    tim.async_wait(std::move(complete));
  }
};

template<typename Timer>
struct test_wait_op_2 : cobalt::op<system::error_code>
{
  Timer & tim;

  test_wait_op_2(Timer & tim) : tim(tim) {}

  void ready(cobalt::handler<system::error_code> h)
  {
    if (tim.expiry() < Timer::clock_type::now())
      h(system::error_code(asio::error::operation_aborted));
  }
  void initiate(cobalt::completion_handler<system::error_code> complete)
  {
    tim.async_wait(std::move(complete));
  }
};


struct post_op : cobalt::op<>
{
  asio::any_io_executor exec;

  post_op(asio::any_io_executor exec) : exec(exec) {}

  void initiate(cobalt::completion_handler<> complete)
  {
    asio::post(std::move(complete));
  }
};


CO_TEST_CASE(op)
{

  asio::steady_timer tim{co_await asio::this_coro::executor, std::chrono::milliseconds(10)};

  co_await test_wait_op{tim};
  co_await test_wait_op{tim};

  tim.expires_after(std::chrono::milliseconds(10));

  co_await test_wait_op_2{tim};
  BOOST_CHECK_THROW(co_await test_wait_op_2{tim}, boost::system::system_error);

  (co_await cobalt::as_result(post_op(co_await asio::this_coro::executor))).value();
  (co_await cobalt::as_result(tim.async_wait(cobalt::use_op))).value();
}

struct op_throw_op
{
  template<typename Handler>
  void operator()(Handler &&)
  {
    throw std::runtime_error("test-exception");

  }
};

template<typename CompletionToken>
auto op_throw(CompletionToken&& token)
{
  return asio::async_initiate<CompletionToken, void(std::exception_ptr)>(
      op_throw_op{}, token);
}


BOOST_AUTO_TEST_CASE(op_throw_)
{

  auto val = [&]() -> cobalt::task<void> {BOOST_CHECK_THROW(co_await op_throw(cobalt::use_op), boost::system::system_error);};

  asio::io_context ctx;
  cobalt::this_thread::set_executor(ctx.get_executor());
  BOOST_CHECK_NO_THROW(spawn(ctx, val(), asio::detached));

BOOST_CHECK_NO_THROW(ctx.run());
}

struct throw_op : cobalt::op<std::exception_ptr>
{
  asio::any_io_executor exec;

  throw_op(asio::any_io_executor exec) : exec(exec) {}

  void initiate(cobalt::completion_handler<std::exception_ptr> complete)
  {
    asio::post(exec, asio::append(std::move(complete), std::make_exception_ptr(std::runtime_error("test-exception"))));
  }
};


CO_TEST_CASE(exception_op)
try
{
  BOOST_CHECK_THROW(co_await throw_op(co_await asio::this_coro::executor), boost::system::system_error);
}
catch(...) {}


struct initiate_op : cobalt::op<>
{
  asio::any_io_executor exec;

  initiate_op(asio::any_io_executor exec) : exec(exec) {}

  void initiate(cobalt::completion_handler<> complete)
  {
    throw std::runtime_error("test-exception");
    asio::post(exec, std::move(complete));
  }
};


CO_TEST_CASE(initiate_exception_op)
try
{
  BOOST_CHECK_THROW(co_await throw_op(co_await asio::this_coro::executor), boost::system::system_error);
} catch(...) {}

CO_TEST_CASE(immediate_executor)
{
  auto called = false;
  asio::post(co_await asio::this_coro::executor, [&]{called = true;});
  asio::experimental::channel<void(system::error_code)> chn{co_await asio::this_coro::executor, 2u};
  co_await chn.async_send(system::error_code(), cobalt::use_op);
  auto [ec] = co_await cobalt::as_tuple(chn.async_receive(cobalt::use_op));
  BOOST_CHECK(!ec);

  BOOST_CHECK(!called);
  co_await cobalt::as_tuple(asio::post(co_await asio::this_coro::executor, cobalt::use_op));
  BOOST_CHECK(called);
}

struct test_async_initiate
{

  template<typename Handler>
  void operator()(Handler && h, std::shared_ptr<int> ptr)
  {
    BOOST_CHECK(ptr);
    asio::dispatch(
        asio::get_associated_immediate_executor(
            h, asio::get_associated_executor(h)),
        std::move(h));
  }
};

template<typename Token>
auto test_cobalt(std::shared_ptr<int> & ptr, Token && token)
{
  return asio::async_initiate<Token, void()>(test_async_initiate{}, token, ptr);
}

CO_TEST_CASE(no_move_from)
{
  std::shared_ptr<int> p = std::make_shared<int>();
  BOOST_CHECK(p);
  co_await test_cobalt(p, cobalt::use_op);
  BOOST_CHECK(p);
}

CO_TEST_CASE(deferred)
{
  asio::steady_timer tim{co_await asio::this_coro::executor, std::chrono::milliseconds(10)};

  co_await post_op(co_await asio::this_coro::executor);
  std::tuple<system::error_code> r = co_await tim.async_wait(asio::as_tuple);
}

template<typename Token>
auto test_multi_values(int a, int b, Token && token)
{
  return asio::async_initiate<Token, void(std::exception_ptr, int, int)>(
      [](auto h, int a, int b)
      {
        asio::dispatch(
            asio::get_associated_immediate_executor(
                h, asio::get_associated_executor(h)),
            [h = std::move(h), a, b]() mutable {
              std::move(h)({}, a, b);
            }
        );
      },
      std::forward<Token>(token),
      a, b
  );
}

CO_TEST_CASE(multi_values)
{
  auto [a, b] = co_await test_multi_values(12, 37, cobalt::use_op);
  BOOST_CHECK_EQUAL(a, 12);
  BOOST_CHECK_EQUAL(b, 37);
}

BOOST_AUTO_TEST_SUITE_END();
