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


static cobalt::promise<std::chrono::milliseconds::rep> dummy(
                                  asio::any_io_executor exec,
                                  std::chrono::milliseconds ms = std::chrono::milliseconds(50))
{
  asio::steady_timer tim{exec, ms};
  co_await tim.async_wait(cobalt::use_op);
  co_return ms.count();
}

static cobalt::generator<int> gen(asio::any_io_executor exec)
{
  asio::steady_timer tim{exec, std::chrono::milliseconds(50000)};
  co_await tim.async_wait(cobalt::use_op);
  co_return 123;
}

BOOST_AUTO_TEST_SUITE(left_race_);
#if !defined(BOOST_COBALT_NO_SELF_DELETE)

CO_TEST_CASE(variadic)
{
  auto exec = co_await asio::this_coro::executor;
  auto d1 = dummy(exec, std::chrono::milliseconds(100));
  auto d2 = dummy(exec, std::chrono::milliseconds( 50));
  auto g = gen(exec);
  auto c = co_await left_race(d1, d2, dummy(exec, std::chrono::milliseconds(100000)), g);
  BOOST_CHECK(c.index() == 1u);
  BOOST_CHECK(boost::variant2::get<1>(c) == 50);
  BOOST_CHECK(d1);
  BOOST_CHECK(!d1.ready());
  BOOST_CHECK( d2.ready());
  BOOST_CHECK(100 == co_await d1);
  BOOST_CHECK(!d1);
  BOOST_CHECK( d1.ready());
  co_await d2;

  g.cancel();
  BOOST_CHECK_THROW(co_await g, boost::system::system_error);
}
#endif

CO_TEST_CASE(list)
{
  auto exec = co_await asio::this_coro::executor;
  std::vector<cobalt::promise<std::chrono::milliseconds::rep>> vec;
  vec.push_back(dummy(exec, std::chrono::milliseconds(100)));
  vec.push_back(dummy(exec, std::chrono::milliseconds( 50)));
  vec.push_back(dummy(exec, std::chrono::milliseconds(100000)));

  auto c = co_await left_race(vec);
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

CO_TEST_CASE(empty_list)
{
  auto exec = co_await asio::this_coro::executor;
  std::vector<cobalt::promise<std::size_t>> vec;
  try {
    BOOST_CHECK_THROW(co_await left_race(vec), boost::system::system_error);
  }
  catch(...) {}
}


CO_TEST_CASE(stop_)
{
  auto d = dummy(co_await asio::this_coro::executor,
                 std::chrono::milliseconds(10));
  BOOST_CHECK((co_await left_race(d, stop())).index() == 0);
}

CO_TEST_CASE(compliance)
{
  auto exec = co_await asio::this_coro::executor;
  auto d = dummy(exec, std::chrono::milliseconds(100000));

  {
    immediate i;
    BOOST_CHECK((co_await left_race(d, i)).index() == 1);
  }

  {
    immediate_bool i;
    BOOST_CHECK((co_await left_race(d, i)).index() == 1);
  }

  {
    immediate_handle i;
    BOOST_CHECK((co_await left_race(d, i)).index() == 1);
  }

  {
    posted p;
    BOOST_CHECK((co_await left_race(d, p)).index() == 1);
  }

  {
    posted_bool p;
    BOOST_CHECK((co_await left_race(d, p)).index() == 1);
  }

  {
    posted_handle p;
    BOOST_CHECK((co_await left_race(d, p)).index() == 1);
  }
  d.cancel();
  BOOST_CHECK_THROW(co_await d, boost::system::system_error);
}
CO_TEST_CASE(immediate_timer)
{
  immediate i;
  asio::steady_timer tim{co_await cobalt::this_coro::executor, std::chrono::steady_clock::time_point::max()};
  BOOST_CHECK((co_await left_race(tim.async_wait(cobalt::use_op), i)) == 1);

}

CO_TEST_CASE(compliance_ranged)
{
  BOOST_CHECK(co_await cobalt::left_race(std::vector<immediate>(3u))        == 0);
  BOOST_CHECK(co_await cobalt::left_race(std::vector<immediate_bool>(1u))   == 0);
  BOOST_CHECK(co_await cobalt::left_race(std::vector<immediate_handle>(1u)) == 0);
  BOOST_CHECK(co_await cobalt::left_race(std::vector<posted>(3u))           == 0);
  BOOST_CHECK(co_await cobalt::left_race(std::vector<posted_bool>(1u))      == 0);
  BOOST_CHECK(co_await cobalt::left_race(std::vector<posted_handle>(1u))    == 0);
}

BOOST_AUTO_TEST_SUITE_END();
