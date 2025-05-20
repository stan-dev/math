//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/race.hpp>
#include <boost/cobalt/generator.hpp>
#include <boost/cobalt/promise.hpp>
#include <boost/cobalt/op.hpp>

#include <boost/asio/steady_timer.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"

using namespace boost;


static  cobalt::promise<std::chrono::milliseconds::rep> dummy(
                                  asio::any_io_executor exec,
                                  std::chrono::milliseconds ms = std::chrono::milliseconds(50))
{
  asio::steady_timer tim{exec, ms};
  co_await tim.async_wait(cobalt::use_op);
  co_return ms.count();
}


static  cobalt::promise<std::chrono::milliseconds::rep> nothrow_dummy(
                                  asio::any_io_executor exec,
                                  std::chrono::milliseconds ms = std::chrono::milliseconds(50))
try {
  asio::steady_timer tim{exec, ms};
  co_await tim.async_wait(cobalt::use_op);
  co_return ms.count();
}
catch(...)
{
  co_return std::numeric_limits<std::chrono::milliseconds::rep>::max();
}

static cobalt::generator<int> gen(asio::any_io_executor exec)
{
  asio::steady_timer tim{exec, std::chrono::milliseconds(50000)};
  co_await tim.async_wait(cobalt::use_op);
  co_return 123;
}

BOOST_AUTO_TEST_SUITE(race_);

CO_TEST_CASE(variadic)
{
  auto exec = co_await asio::this_coro::executor;
  auto d1 = dummy(exec, std::chrono::milliseconds(100));
  auto d2 = dummy(exec, std::chrono::milliseconds( 50));
  auto g = gen(exec);
  std::mt19937 src{1u};
  auto c = co_await race(src, d1, d2, dummy(exec, std::chrono::milliseconds(100000)), g);
  BOOST_CHECK(c.index() == 1u);
  BOOST_CHECK(boost::variant2::get<1>(c) == 50);
  BOOST_CHECK(d1);
  //BOOST_CHECK(!d1.ready()); NOTE: Inderministic on msvc, due to the additional post!
  BOOST_CHECK( d2.ready());
  BOOST_CHECK(100 == co_await d1);
  BOOST_CHECK(!d1);
  BOOST_CHECK( d1.ready());
  co_await d2;

  g.cancel();
  try { BOOST_CHECK_THROW(co_await g, boost::system::system_error); }
  catch (boost::system::system_error &) {}
}


cobalt::promise<void> list_step(std::default_random_engine::result_type seed)
{
  std::mt19937 src{seed};

  auto exec = co_await asio::this_coro::executor;
  std::vector<cobalt::promise<std::chrono::milliseconds::rep>> vec;
  vec.push_back(dummy(exec, std::chrono::milliseconds(100)));
  vec.push_back(dummy(exec, std::chrono::milliseconds( 50)));
  vec.push_back(dummy(exec, std::chrono::milliseconds(100000)));

  auto c = co_await race(src, vec);
  BOOST_CHECK(c.first == 1u);
  BOOST_CHECK(c.second == 50);
  BOOST_CHECK(!vec[0].ready());
  BOOST_CHECK( vec[1].ready());
  BOOST_CHECK(co_await vec[0]);
  BOOST_CHECK( vec[0].ready());
  BOOST_CHECK( vec[1].ready());
  vec[2].cancel();
  BOOST_CHECK( vec[2]);
  BOOST_CHECK_THROW(co_await vec[2], boost::system::system_error);
  BOOST_CHECK_THROW(co_await vec[2], boost::system::system_error);
  BOOST_CHECK(!vec[2]);
}

CO_TEST_CASE(list_1u) { co_await list_step(1u);}
CO_TEST_CASE(list_2u) { co_await list_step(2u);}
CO_TEST_CASE(list_3u) { co_await list_step(3u);}
CO_TEST_CASE(list_4u) { co_await list_step(4u);}
CO_TEST_CASE(list_5u) { co_await list_step(5u);}
CO_TEST_CASE(list_6u) { co_await list_step(6u);}
CO_TEST_CASE(list_7u) { co_await list_step(7u);}
CO_TEST_CASE(list_8u) { co_await list_step(8u);}
CO_TEST_CASE(list_9u) { co_await list_step(9u);}

CO_TEST_CASE(empty_list)
{
  auto exec = co_await asio::this_coro::executor;
  std::vector<cobalt::promise<std::size_t>> vec;
  try
  {
    BOOST_CHECK_THROW(co_await race(vec),  boost::system::system_error);
  }
  catch(...) {}
}


CO_TEST_CASE(stop_)
{
  auto d = nothrow_dummy(co_await asio::this_coro::executor,
                 std::chrono::milliseconds(10));
  BOOST_CHECK((co_await left_race(d, stop())).index() == 0);
}

CO_TEST_CASE(compliance)
{
  auto exec = co_await asio::this_coro::executor;
  auto d = dummy(exec, std::chrono::milliseconds(100000));
  {
    immediate i;
    BOOST_CHECK((co_await race(d, i)).index() == 1);
  }

  {
    immediate_bool i;
    BOOST_CHECK((co_await race(d, i)).index() == 1);
  }

  {
    immediate_handle i;
    BOOST_CHECK((co_await race(d, i)).index() == 1);
  }
  {
    posted p;
    BOOST_CHECK((co_await race(d, p)).index() == 1);
  }
  {
    posted_bool p;
    BOOST_CHECK((co_await race(d, p)).index() == 1);
  }
  {
    posted_handle p;
    BOOST_CHECK((co_await race(d, p)).index() == 1);
  }
  d.cancel();
  BOOST_CHECK_THROW(co_await d, boost::system::system_error);
}


BOOST_AUTO_TEST_SUITE_END();
