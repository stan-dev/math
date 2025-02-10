//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/thread.hpp>

#include <boost/asio/append.hpp>
#include <boost/asio/bind_executor.hpp>

#include <boost/core/no_exceptions_support.hpp>

namespace boost::cobalt
{

namespace detail
{

thread_promise::thread_promise()
    : promise_cancellation_base<asio::cancellation_slot, asio::enable_total_cancellation>(
    signal_helper_2::signal.slot(), asio::enable_total_cancellation())
{
  mtx.lock();
}

void run_thread(
    std::shared_ptr<thread_state> st_,
    unique_handle<thread_promise> h)
{

#if !defined(BOOST_COBALT_NO_PMR)
  pmr::unsynchronized_pool_resource resource;
  boost::cobalt::this_thread::set_default_resource(&resource);
  h->resource = &resource;
#endif

  {
    auto st = std::move(st_);
    h->reset_cancellation_source(st->signal.slot());
    h->set_executor(st->ctx.get_executor());
    boost::cobalt::this_thread::set_executor(st->ctx.get_executor());

    asio::post(
        st->ctx.get_executor(),
        [st, h = std::move(h)]() mutable
        {
          std::lock_guard<std::mutex> lock{h->mtx};
          std::move(h).resume();
        });

    std::exception_ptr ep;

    BOOST_TRY
    {
      st->ctx.run();
    }
    BOOST_CATCH(...)
    {
      ep = std::current_exception();
    }
    BOOST_CATCH_END

    st->done = true;
    st->signal.slot().clear();
    std::lock_guard<std::mutex> lock(st->mtx);
    if (!st->waitor && ep) // nobodies waiting, so unhandled exception
      std::rethrow_exception(ep);
    else if (st->waitor)
      asio::post(asio::append(*std::exchange(st->waitor, std::nullopt), ep));
  }
}


boost::cobalt::thread detail::thread_promise::get_return_object()
{
  auto st = std::make_shared<thread_state>();
  boost::cobalt::thread res{std::thread{
      [st,
       h = unique_handle<detail::thread_promise>::from_promise(*this)]() mutable
      {
        run_thread(std::move(st), std::move(h));
      }
     }, st
    };

  return res;
}


}

void thread::join() {thread_.join();}
bool thread::joinable() const {return thread_.joinable();}
void thread::detach()
{
  thread_.detach();
  state_ = nullptr;
}
}