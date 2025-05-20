// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/strand.hpp>
#include <boost/asio/io_context_strand.hpp>

#if defined(BOOST_COBALT_USE_BOOST_CONTAINER_PMR)
#include <boost/container/pmr/synchronized_pool_resource.hpp>
#endif

#include <boost/test/unit_test.hpp>
#include "test.hpp"

using namespace boost;


BOOST_AUTO_TEST_SUITE(strand);

cobalt::promise<void> do_the_thing()
{
  static std::atomic<int> unique = 0;

  asio::steady_timer tim{co_await cobalt::this_coro::executor, std::chrono::milliseconds(50)};
  co_await tim.async_wait(cobalt::use_op);
  unique++;
  BOOST_CHECK(unique == 1);
  unique--;
}

#if !defined(BOOST_COBALT_USE_IO_CONTEXT)

BOOST_AUTO_TEST_CASE(strand)
{
  std::vector<std::thread> ths;

  asio::io_context ctx;
  asio::any_io_executor exec{asio::make_strand(ctx.get_executor())};
#if !defined(BOOST_COBALT_NO_PMR)
  cobalt::pmr::synchronized_pool_resource sync;
#endif

  for (int i = 0; i < 8; i++)
    ths.push_back(
        std::thread{
          [&]
          {
#if !defined(BOOST_COBALT_NO_PMR)
            cobalt::this_thread::set_default_resource(&sync);
#endif
            cobalt::this_thread::set_executor(exec);
            +do_the_thing();
            ctx.run();
          }});


  for (auto & th : ths)
    th.join();
}

#else

BOOST_AUTO_TEST_CASE(dummy) {}

#endif

BOOST_AUTO_TEST_SUITE_END();

