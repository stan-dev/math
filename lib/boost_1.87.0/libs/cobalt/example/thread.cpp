// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/// This example shows how to use threads to offload cpu_intense work.

#include <boost/cobalt.hpp>

#include <boost/asio/as_tuple.hpp>
#include <boost/asio/redirect_error.hpp>
#include <boost/asio/experimental/concurrent_channel.hpp>
#include <boost/asio/this_coro.hpp>


namespace cobalt = boost::cobalt;
using boost::system::error_code;

template<typename Signature>
using cchannel = boost::asio::experimental::concurrent_channel<Signature>;

// this is a function doing some CPU heavy work that should be offloaded onto a thread
cobalt::promise<int> cpu_intense_work(int a, int b) {co_return a + b;}

// this channel is used to send a response to completed work
using response_channel = cchannel<void(std::exception_ptr, int)>;
// this channel is used to send a request to a working thread
using request_channel = cchannel<void(error_code, int, int, response_channel * res)>;

// the worker wrapper
cobalt::thread worker(request_channel & work)
{
  while (work.is_open())
  {
    auto [ec, a, b, respond_to] = co_await work.async_receive(boost::asio::as_tuple(cobalt::use_op));
    if (ec) // done, ignore. in our code this is only triggered by closing the channel
      break;

    // to emulate this being like awaiting on the same thread, we also deliver an exception.
    std::exception_ptr ep;
    int res = 0;
    try
    {
      res = co_await cpu_intense_work(a, b);
    }
    catch(...)
    {
      // this way exception get sent to the awaiting coro as if it was a call.
      ep = std::current_exception();
    }
    // send the response. If the channel is closed, the program will terminate!
    co_await respond_to->async_send(ep, res, boost::asio::redirect_error(cobalt::use_op, ec));
  }
}

cobalt::promise<void> work(request_channel & rc, int min_a, int max_a, int b)
{
  response_channel res{co_await cobalt::this_coro::executor};
  for (int a = min_a; a <= max_a; a++)
  {
    // the following two calls offload the work to another thread.
    co_await rc.async_send(error_code{}, a, b, &res, cobalt::use_op);
    int c = co_await res.async_receive(cobalt::use_op); // may throw if working thread has an exception
    printf("The CPU intensive result of adding %d to %d, is %d\n", a, b, c);
  }
}

cobalt::main co_main(int argc, char *argv [])
{
  // a very simple thread pool
  std::vector<cobalt::thread> thrs;
  const std::size_t n = 4u;

  request_channel rc{co_await cobalt::this_coro::executor};
  for (auto i = 0u; i < n; i++)
      thrs.push_back(worker(rc));

  try
  {
    // this is an over simplification, but emulated multiple pieces of
    // code in the single threaded environment offloading work to the thread.
    co_await cobalt::join(
        work(rc, 0, 10, 32),
        work(rc, 10, 20, 22),
        work(rc, 50, 60, -18)
        );

  }
  catch(std::exception & e)
  {
    printf("Completed with exception %s\n", e.what());
  }
  // closing the channel will cause the threads to complete
  rc.close();
  // wait them so they don't leak.
  co_await cobalt::join(thrs);
  co_return 0;
}