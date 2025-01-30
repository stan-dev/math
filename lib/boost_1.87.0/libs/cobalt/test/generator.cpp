//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/generator.hpp>
#include <boost/cobalt/promise.hpp>
#include <boost/cobalt/race.hpp>
#include <boost/cobalt/op.hpp>
#include <boost/core/ignore_unused.hpp>

#include <boost/asio/steady_timer.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"

using namespace boost;

BOOST_AUTO_TEST_SUITE(generator);

cobalt::generator<int> gen()
{
  for (int i = 0; i <10; i ++)
    co_yield i;

  co_return 10;
}



CO_TEST_CASE(generator_int)
{
  auto g = gen();
  int i = 0;
  {
    auto aw = g.operator co_await();
    BOOST_CHECK(aw.await_ready());
    BOOST_CHECK(i ++ == co_await std::move(aw));
  }

  while (g)
    BOOST_CHECK(i ++ == co_await g);



  BOOST_CHECK(i == 11);

  g = gen();
  BOOST_CHECK(g);
  while (g)
    co_await g;

  co_return ;
}


cobalt::generator<int, int> gen_push()
{
  int val = 1u;
  for (int i = 0; i < 10; i++)
  {
    auto v = co_yield val;
    BOOST_CHECK(v == val);
    val += v;
  }


  co_return val;
}


CO_TEST_CASE(generator_push)
{
  auto g = gen_push();

  int i = 1;
  int nw = 1;
  while (g)
  {
    nw = co_await g(i);
    BOOST_CHECK(i == nw);
    i *= 2;
  }

  BOOST_CHECK(i == 2048);
  co_return ;
}

cobalt::generator<int> delay_gen(std::chrono::milliseconds tick)
{
  asio::steady_timer tim{co_await cobalt::this_coro::executor, std::chrono::steady_clock::now()};
  for (int i = 0; i < 10; i ++)
  {
    co_await tim.async_wait(cobalt::use_op);
    tim.expires_at(tim.expiry() + tick);
    co_yield i;
  }
  co_return 10;
}

cobalt::generator<int> lazy_delay_gen(std::chrono::milliseconds tick)
{
  co_await cobalt::this_coro::initial;
  asio::steady_timer tim{co_await cobalt::this_coro::executor, std::chrono::steady_clock::now()};
  for (int i = 0; i < 10; i ++)
  {
    co_await tim.async_wait(cobalt::use_op);
    tim.expires_at(tim.expiry() + tick);
    co_yield i;
  }
  co_return 10;
}


#if !defined(BOOST_COBALT_NO_SELF_DELETE)

CO_TEST_CASE(generator_left_race)
{
  asio::steady_timer tim{co_await cobalt::this_coro::executor, std::chrono::milliseconds(50)};
  auto g1 = delay_gen(std::chrono::milliseconds(200));
  co_await tim.async_wait(cobalt::use_op);
  auto g2 = delay_gen(std::chrono::milliseconds(100));

  using v = variant2::variant<int, int>;
  auto v1 = [](int value) -> v {return v{variant2::in_place_index<0u>, value};};
  auto v2 = [](int value) -> v {return v{variant2::in_place_index<1u>, value};};

  BOOST_CHECK(v1(0) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(0) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(1) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(1) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(2) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(3) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(2) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(4) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(5) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(3) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(6) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(7) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(4) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(8) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(9) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(10) == co_await left_race(g1, g2));


  BOOST_CHECK(!g2);
  g1.cancel();
  BOOST_CHECK_THROW(co_await g1, boost::system::system_error);
}


CO_TEST_CASE(lazy_generator_left_race)
{
  asio::steady_timer tim{co_await cobalt::this_coro::executor, std::chrono::milliseconds(50)};
  auto g1 = lazy_delay_gen(std::chrono::milliseconds(200));
  co_await tim.async_wait(cobalt::use_op);
  auto g2 = lazy_delay_gen(std::chrono::milliseconds(100));

  using v = variant2::variant<int, int>;
  auto v1 = [](int value) -> v {return v{variant2::in_place_index<0u>, value};};
  auto v2 = [](int value) -> v {return v{variant2::in_place_index<1u>, value};};

  BOOST_CHECK(v1(0) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(0) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(1) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(2) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(1) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(3) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(4) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(2) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(5) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(6) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(3) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(7) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(8) == co_await left_race(g1, g2));
  BOOST_CHECK(v1(4) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(9) == co_await left_race(g1, g2));
  BOOST_CHECK(v2(10) == co_await left_race(g1, g2));


  BOOST_CHECK(!g2);
  g1.cancel();
  BOOST_CHECK_THROW(co_await g1, boost::system::system_error);
}

#endif

cobalt::generator<int> gshould_unwind(asio::io_context & ctx)
{
  co_await asio::post(ctx, cobalt::use_op);
  co_return 0;
}

BOOST_AUTO_TEST_CASE(unwind)
{
  asio::io_context ctx;
  boost::cobalt::this_thread::set_executor(ctx.get_executor());
  boost::ignore_unused(gshould_unwind(ctx));
}


cobalt::generator<int> gen_stop()
{
  int val = 1u;
  for (int i = 0; i < 10; i++)
  {
    if (i == 4)
      co_await stop();
    co_yield i;

  }
  co_return val;
}

#if !defined(BOOST_COBALT_USE_BOOST_CONTAINER_PMR)
// clang-14 does not like this test.

CO_TEST_CASE(stop)
{
  auto g = gen_stop();
  while (g)
    co_await g;

  auto gg = std::move(g);
}

#endif

cobalt::generator<int, int> eager()
{
  int i = co_await cobalt::this_coro::initial;
  for (; i < 10; i += co_yield i);

  co_return i;
}


CO_TEST_CASE(eager_2)
{
  auto g = eager();


  BOOST_CHECK(2 == co_await g(2));
  BOOST_CHECK(6 == co_await g(4));
  BOOST_CHECK(9 == co_await g(3));

  auto gg = std::move(g);
//  BOOST_CHECK(!g);
  BOOST_CHECK(gg);
  BOOST_CHECK(15 == co_await gg(6));
  BOOST_CHECK(!gg);

}

struct generator_move_only
{
  generator_move_only() = default;
  generator_move_only(generator_move_only &&) = default;
  generator_move_only & operator=(generator_move_only &&) = delete;
};

cobalt::generator<generator_move_only, generator_move_only> gen_move_only_test()
{
  co_yield generator_move_only{};
  co_return generator_move_only{};
}

CO_TEST_CASE(move_only)
{
  auto g = gen_move_only_test();

  co_await g(generator_move_only{});
  co_await g(generator_move_only{});
}

cobalt::generator<int> detached()
{
  co_await asio::post(cobalt::use_op);
  co_yield 42;
  co_return 24;
}

CO_TEST_CASE(detached_)
{
  boost::ignore_unused(detached());
  co_return;
}

cobalt::generator<int, int> detached_push()
{
  int i = co_yield 42;
  while (true)
    i = co_yield i;
  co_return i;
}

cobalt::generator<int> np() {return cobalt::noop(42);}

CO_TEST_CASE(detached_push_)
{
  auto g = detached_push();

  co_await g(1);
  co_await np();
}

BOOST_AUTO_TEST_SUITE_END();