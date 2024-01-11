//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/gather.hpp>
#include <boost/cobalt/generator.hpp>
#include <boost/cobalt/promise.hpp>
#include <boost/cobalt/op.hpp>

#include <boost/asio/steady_timer.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"

using namespace boost;

static cobalt::promise<std::chrono::milliseconds::rep> wdummy(
                                          asio::any_io_executor exec,
                                          std::chrono::milliseconds ms = std::chrono::milliseconds(25))
{
  if (ms == std::chrono::milliseconds ::max())
    throw std::runtime_error("wdummy_throw");

  asio::steady_timer tim{exec, ms};
  co_await tim.async_wait(cobalt::use_op);
  co_return ms.count();
}

static cobalt::generator<int> wgen(asio::any_io_executor exec)
{
  asio::steady_timer tim{exec, std::chrono::milliseconds(25)};
  co_await tim.async_wait(cobalt::use_op);
  co_return 123;
}

static cobalt::promise<void> wthrow()
{
  throw std::runtime_error("wthrow");
  co_return;
}


BOOST_AUTO_TEST_SUITE(gather_);

CO_TEST_CASE(variadic)
{
  auto exec = co_await asio::this_coro::executor;
  auto d1 = wdummy(exec, std::chrono::milliseconds(100));
  auto d2 = wdummy(exec, std::chrono::milliseconds( 50));
  asio::steady_timer tim{co_await asio::this_coro::executor};
  auto g = wgen(exec);
  auto c = co_await cobalt::gather(d1, d2, wdummy(exec, std::chrono::milliseconds(150)),
                           g, wthrow());

  BOOST_CHECK( std::get<0>(c).has_value());
  BOOST_CHECK(!std::get<0>(c).has_error());
  BOOST_CHECK(*std::get<0>(c) == 100);
  BOOST_CHECK( std::get<2>(c).has_value());
  BOOST_CHECK(!std::get<1>(c).has_error());
  BOOST_CHECK(*std::get<1>(c) ==  50);
  BOOST_CHECK( std::get<2>(c).has_value());
  BOOST_CHECK(!std::get<2>(c).has_error());
  BOOST_CHECK(*std::get<2>(c) ==  150);
  BOOST_CHECK( std::get<3>(c).has_value());
  BOOST_CHECK(!std::get<3>(c).has_error());
  BOOST_CHECK(*std::get<3>(c) ==  123);
  BOOST_CHECK(!std::get<4>(c).has_value());
  BOOST_CHECK( std::get<4>(c).has_error());
  try {
  BOOST_CHECK_THROW(std::get<4>(c).value(), boost::system::system_error);
  } catch(...) {}
}

CO_TEST_CASE(list)
{
  auto exec = co_await asio::this_coro::executor;
  std::vector<cobalt::promise<std::chrono::milliseconds::rep>> vec;
  vec.push_back(wdummy(exec, std::chrono::milliseconds(100)));
  vec.push_back(wdummy(exec, std::chrono::milliseconds( 50)));
  vec.push_back(wdummy(exec, std::chrono::milliseconds(150)));
  vec.push_back(wdummy(exec, std::chrono::milliseconds::max()));

  auto res = co_await gather(vec);
  BOOST_REQUIRE(res.size() == 4);
  BOOST_CHECK( res[0].has_value());
  BOOST_CHECK(!res[0].has_error());
  BOOST_CHECK( res[0].value() == 100);
  BOOST_CHECK( res[1].has_value());
  BOOST_CHECK(!res[1].has_error());
  BOOST_CHECK( res[1].value() == 50);
  BOOST_CHECK( res[2].has_value());
  BOOST_CHECK(!res[2].has_error());
  BOOST_CHECK( res[2].value() == 150);
  BOOST_CHECK(!res[3].has_value());
  BOOST_CHECK( res[3].has_error());
  try {
    BOOST_CHECK_THROW(res[3].value(), boost::system::system_error);
  } catch(...) {}
}

CO_TEST_CASE(compliance)
{
  auto exec = co_await asio::this_coro::executor;
  {
    auto d = wdummy(exec, std::chrono::milliseconds(1));
    immediate i;
    co_await gather(d, i);
  }

  {
    auto d = wdummy(exec, std::chrono::milliseconds(1));
    immediate_bool i;
    co_await gather(d, i);
  }
  {
    auto d = wdummy(exec, std::chrono::milliseconds(1));
    immediate_handle i;
    co_await gather(d, i);
  }
  {
    auto d = wdummy(exec, std::chrono::milliseconds(1));
    posted p;
    co_await gather(d, p);
  }
  {
    auto d = wdummy(exec, std::chrono::milliseconds(1));
    posted_bool p;
    co_await gather(d, p);
  }
  {
    auto d = wdummy(exec, std::chrono::milliseconds(1));
    posted_handle p;
    co_await gather(d, p);
  }
}

BOOST_AUTO_TEST_SUITE_END();