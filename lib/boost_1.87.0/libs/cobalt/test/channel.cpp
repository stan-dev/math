//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/channel.hpp>
#include <boost/cobalt/promise.hpp>
#include <boost/cobalt/race.hpp>
#include <boost/cobalt/gather.hpp>
#include <boost/cobalt/async_for.hpp>

#include <boost/cobalt/join.hpp>
#include <boost/asio/steady_timer.hpp>

#include "test.hpp"
#include <boost/test/unit_test.hpp>
#include <any>

namespace cobalt = boost::cobalt;

cobalt::promise<void> do_write(cobalt::channel<void> &chn, std::vector<int> & seq)
{
  seq.push_back(0);
  co_await chn.write(); seq.push_back(1);
  co_await chn.write(); seq.push_back(2);
  (co_await cobalt::as_result(chn.write())).value(); seq.push_back(3);
  co_await cobalt::as_tuple(chn.write()); seq.push_back(4);
  co_await chn.write(); seq.push_back(5);
  co_await chn.write(); seq.push_back(6);
  co_await chn.write(); seq.push_back(7);
}

cobalt::promise<void> do_read(cobalt::channel<void> &chn, std::vector<int> & seq)
{
  seq.push_back(10);
  co_await chn.read(); seq.push_back(11);
  co_await chn.read(); seq.push_back(12);
  (co_await cobalt::as_result(chn.read())).value(); seq.push_back(13);
  co_await cobalt::as_tuple(chn.read()); seq.push_back(14);
  co_await chn.read(); seq.push_back(15);
  co_await chn.read(); seq.push_back(16);
  co_await chn.read(); seq.push_back(17);
}

BOOST_AUTO_TEST_SUITE(channel);

CO_TEST_CASE(void_)
{
  cobalt::channel<void> chn{2u, co_await cobalt::this_coro::executor};

  std::vector<int> seq;
  auto r = do_read(chn, seq);
  auto w = do_write(chn, seq);

  co_await r;
  co_await w;
  BOOST_REQUIRE(seq.size() == 16);
  BOOST_CHECK(seq[0] == 10);
  BOOST_CHECK(seq[1] == 0);
  BOOST_CHECK(seq[2] == 1);
  BOOST_CHECK(seq[3] == 2);
  BOOST_CHECK(seq[4] == 11);
  BOOST_CHECK(seq[5] == 12);
  BOOST_CHECK(seq[6] == 3);
  BOOST_CHECK(seq[7] == 4);
  BOOST_CHECK(seq[8] == 13);
  BOOST_CHECK(seq[9] == 14);
  BOOST_CHECK(seq[10] == 5);
  BOOST_CHECK(seq[11] == 6);
  BOOST_CHECK(seq[12] == 15);
  BOOST_CHECK(seq[13] == 16);
  BOOST_CHECK(seq[14] == 7);
  BOOST_CHECK(seq[15] == 17);
}

CO_TEST_CASE(void_0)
{
    cobalt::channel<void> chn{0u, co_await cobalt::this_coro::executor};

    std::vector<int> seq;
    auto r = do_read(chn, seq);
    auto w = do_write(chn, seq);

    co_await r;
    co_await w;
    BOOST_REQUIRE(seq.size() == 16);
    BOOST_CHECK(seq[0] == 10);
    BOOST_CHECK(seq[1] == 0);
    BOOST_CHECK(seq[2] == 11);
    BOOST_CHECK(seq[3] == 1);
    BOOST_CHECK(seq[4] == 12);
    BOOST_CHECK(seq[5] == 2);
    BOOST_CHECK(seq[6] == 13);
    BOOST_CHECK(seq[7] == 3);
    BOOST_CHECK(seq[8] == 14);
    BOOST_CHECK(seq[9] == 4);
    BOOST_CHECK(seq[10] == 15);
    BOOST_CHECK(seq[11] == 5);
    BOOST_CHECK(seq[12] == 16);
    BOOST_CHECK(seq[13] == 6);
    BOOST_CHECK(seq[14] == 17);
    BOOST_CHECK(seq[15] == 7);
}

cobalt::promise<void> do_write(cobalt::channel<int> &chn, std::vector<int> & seq)
{
  seq.push_back(0);
  co_await chn.write(1); seq.push_back(1);
  co_await chn.write(2); seq.push_back(2);
  (co_await cobalt::as_result(chn.write(3))).value(); seq.push_back(3);
  co_await cobalt::as_tuple(chn.write(4)); seq.push_back(4);
  co_await chn.write(5); seq.push_back(5);
  co_await chn.write(6); seq.push_back(6);
  co_await chn.write(7); seq.push_back(7);
}

cobalt::promise<void> do_read(cobalt::channel<int> &chn, std::vector<int> & seq)
{
  seq.push_back(10);
  BOOST_CHECK(1 == co_await chn.read()); seq.push_back(11);
  BOOST_CHECK(2 == co_await chn.read()); seq.push_back(12);
  BOOST_CHECK(3 == (co_await cobalt::as_result(chn.read())).value()); seq.push_back(13);
  BOOST_CHECK(4 == std::get<1>(co_await cobalt::as_tuple(chn.read()))); seq.push_back(14);
  BOOST_CHECK(5 == co_await chn.read()); seq.push_back(15);
  BOOST_CHECK(6 == co_await chn.read()); seq.push_back(16);
  BOOST_CHECK(7 == co_await chn.read()); seq.push_back(17);
}


CO_TEST_CASE(int_)
{
  cobalt::channel<int> chn{2u, co_await cobalt::this_coro::executor};

  std::vector<int> seq;
  auto w = do_write(chn, seq);
  auto r = do_read(chn, seq);

  co_await r;
  co_await w;
  BOOST_REQUIRE(seq.size() == 16);
  BOOST_CHECK(seq[0] == 0);
  BOOST_CHECK(seq[1] == 1);
  BOOST_CHECK(seq[2] == 2);
  BOOST_CHECK(seq[3] == 10);
  BOOST_CHECK(seq[4] == 11);
  BOOST_CHECK(seq[5] == 12);
  BOOST_CHECK(seq[6] == 3);
  BOOST_CHECK(seq[7] == 4);
  BOOST_CHECK(seq[8] == 13);
  BOOST_CHECK(seq[9] == 14);
  BOOST_CHECK(seq[10] == 5);
  BOOST_CHECK(seq[11] == 6);
  BOOST_CHECK(seq[12] == 15);
  BOOST_CHECK(seq[13] == 16);
  BOOST_CHECK(seq[14] == 7);
  BOOST_CHECK(seq[15] == 17);
}

cobalt::promise<void> do_write(cobalt::channel<std::string> &chn, std::vector<int> & seq)
{
    seq.push_back(0);
    co_await chn.write(std::string("1")); seq.push_back(1);
    co_await chn.write(std::string("2")); seq.push_back(2);
    co_await chn.write(std::string("3")); seq.push_back(3);
    co_await chn.write(std::string("4")); seq.push_back(4);
    co_await chn.write(std::string("5")); seq.push_back(5);
    co_await chn.write(std::string("6")); seq.push_back(6);
    co_await chn.write(std::string("7 but we need to be sure we get ouf of SSO")); seq.push_back(7);
}

cobalt::promise<void> do_read(cobalt::channel<std::string> &chn, std::vector<int> & seq)
{
    seq.push_back(10);
    BOOST_CHECK("1" == co_await chn.read()); seq.push_back(11);
    BOOST_CHECK("2" == co_await chn.read()); seq.push_back(12);
    BOOST_CHECK("3" == co_await chn.read()); seq.push_back(13);
    BOOST_CHECK("4" == co_await chn.read()); seq.push_back(14);
    BOOST_CHECK("5" == co_await chn.read()); seq.push_back(15);
    BOOST_CHECK("6" == co_await chn.read()); seq.push_back(16);
    BOOST_CHECK("7 but we need to be sure we get ouf of SSO" == co_await chn.read()); seq.push_back(17);
}


CO_TEST_CASE(str)
{
    cobalt::channel<std::string> chn{0u, co_await cobalt::this_coro::executor};

    std::vector<int> seq;
    auto w = do_write(chn, seq);
    auto r = do_read(chn, seq);

    co_await r;
    co_await w;
    BOOST_REQUIRE(seq.size() == 16);
    BOOST_CHECK(seq[0] == 0);
    BOOST_CHECK(seq[1] == 10);
    BOOST_CHECK(seq[2] == 1);
    BOOST_CHECK(seq[3] == 11);
    BOOST_CHECK(seq[4] == 2);
    BOOST_CHECK(seq[5] == 12);
    BOOST_CHECK(seq[6] == 3);
    BOOST_CHECK(seq[7] == 13);
    BOOST_CHECK(seq[8] == 4);
    BOOST_CHECK(seq[9] == 14);
    BOOST_CHECK(seq[10] == 5);
    BOOST_CHECK(seq[11] == 15);
    BOOST_CHECK(seq[12] == 6);
    BOOST_CHECK(seq[13] == 16);
    BOOST_CHECK(seq[14] == 7);
    BOOST_CHECK(seq[15] == 17);
}

CO_TEST_CASE(raceable)
{
    cobalt::channel<int>  ci{0u};
    cobalt::channel<void> cv{0u};
    auto [r1, r2] = co_await cobalt::gather(cobalt::race(ci.read(), cv.read()), cv.write());
    r1.value();
    BOOST_REQUIRE(r1.has_value());
    BOOST_CHECK(r1->index() == 1u);
    BOOST_CHECK(!r2.has_error());
}

CO_TEST_CASE(raceable_1)
{
  cobalt::channel<int>  ci{1u};
  cobalt::channel<void> cv{1u};
  auto [r1, r2] = co_await cobalt::gather(
      cobalt::race(ci.read(), cv.read()),
      cv.write());
  BOOST_CHECK(r1->index() == 1u);
  BOOST_CHECK(!r2.has_error());
}



namespace issue_53
{

cobalt::promise<void> timeout_and_write(cobalt::channel<std::string> &channel)
{
  while (!co_await cobalt::this_coro::cancelled)
  {
    boost::asio::steady_timer timer{co_await cobalt::this_coro::executor};
    timer.expires_after(std::chrono::seconds{20});
    co_await timer.async_wait(cobalt::use_op);
    std::string val("Test!");
    co_await channel.write(val);
  }

  co_return;
}

cobalt::promise<void> read(cobalt::channel<std::string> &channel)
{
  while (!co_await cobalt::this_coro::cancelled)
    co_await channel.read();
}

cobalt::promise<void> test()
{
  cobalt::channel<std::string> channel;
  co_await cobalt::join(timeout_and_write(channel), read(channel));
}

CO_TEST_CASE(issue_93)
{
  co_await cobalt::race(test(), boost::asio::post(cobalt::use_op));
}

cobalt::promise<void> writer(cobalt::channel<int> & c)
{
  for (int i = 0; i < 10; i++)
    co_await c.write(i);
  c.close();
}

CO_TEST_CASE(reader)
{
  cobalt::channel<int> c;

  +writer(c);
  int i = 0;
  BOOST_COBALT_FOR(int value, cobalt::channel_reader(c))
    BOOST_CHECK(value == i++);

}

}


BOOST_AUTO_TEST_SUITE_END();

namespace boost::cobalt
{

CO_TEST_CASE(unique)
{
  std::unique_ptr<int> p{new int(42)};
  auto pi = p.get();
  cobalt::channel<std::unique_ptr<int>> c{1u};

  co_await c.write(std::move(p));
  auto p2 = co_await c.read();

  BOOST_CHECK(p == nullptr);
  BOOST_CHECK(p2.get() == pi);
}

CO_TEST_CASE(any)
{
  cobalt::channel<std::any> c{1u};
  co_return ;
}



}
