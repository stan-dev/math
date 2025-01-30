// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/main.hpp>
#include <boost/cobalt/op.hpp>
#include <boost/asio/steady_timer.hpp>

using namespace boost;

// tag::timer_example[]
struct wait_op final : cobalt::op<system::error_code> // <1>
{
  asio::steady_timer & tim;
  wait_op(asio::steady_timer & tim) : tim(tim) {}
  void ready(cobalt::handler<system::error_code> h ) override // <2>
  {
    if (tim.expiry() < std::chrono::steady_clock::now())
      h(system::error_code{});
  }
  void initiate(cobalt::completion_handler<system::error_code> complete) override // <3>
  {
    tim.async_wait(std::move(complete));
  }
};


cobalt::main co_main(int argc, char * argv[])
{
  asio::steady_timer tim{co_await asio::this_coro::executor,
                         std::chrono::milliseconds(std::stoi(argv[1]))};
  co_await wait_op(tim); // <4>
  co_return 0; //
}
// end::timer_example[]
