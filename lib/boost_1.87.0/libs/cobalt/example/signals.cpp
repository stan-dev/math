//
// Copyright (c) 2023 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//


#include <boost/signals2.hpp>
#include <boost/cobalt.hpp>
#include <boost/callable_traits/args.hpp>

namespace cobalt = boost::cobalt;
namespace signals = boost::signals2;

template<typename Signal>
struct signal_awaitable
{
  using args_type = boost::callable_traits::args_t<typename Signal::signature_type>;


  bool await_ready() { return false; } // < always wait for the signal to fire.
  void await_suspend(std::coroutine_handle<void> h)
  {
    awaited_from.reset(h.address());
    // the handler will get copied, so we can't capture the handle with a unique_ptr
    signal.connect_extended(
        [this, _ = boost::intrusive_ptr<signal_awaitable>(this)
         ](const signals::connection & conn, auto ... args) mutable
        {
          auto aw = std::move(awaited_from);
          conn.disconnect();
          result_cache.emplace(std::move(args)...); // the result_catch lives in the coro frame
          std::move(aw).resume(); // release ownership & resume
        });

  }

  auto await_resume() // return the value
  {
    constexpr std::size_t size = std::tuple_size_v<args_type>;
    if constexpr (size == 1u) // single argument doesn't need a tuple
      return std::get<0u>(*std::move(result_cache));
    else if constexpr (size > 1u) // make a tuple if more than one arg
      return *std::move(result_cache);
    // else return void.
  }

  // capture it for lazy initialization
  Signal & signal;
  // capture the ownership of the awaiting coroutine
  cobalt::unique_handle<void> awaited_from;
  // store the result from the call
  std::optional<args_type> result_cache;

  // to manage shared ownership with an internal counter.
  // If the last gets released before the handler is invoked,
  // the coro just gets destroyed.
  std::size_t use_count{0u};
  friend void intrusive_ptr_add_ref(signal_awaitable * st) {st->use_count++;}
  friend void intrusive_ptr_release(signal_awaitable * st)
  {
    if (st->use_count-- == 1u)
      st->awaited_from.reset();
  }
};

namespace boost::signals2
{

// make all signals awaitable
template<typename ... Args>
auto operator co_await(signals::signal<Args...> & sig) -> signal_awaitable<signals::signal<Args...>>
{
  return {sig};
}

}

cobalt::promise<int> await_signal(signals::signal<void(int)>  & sig)
{
  co_return co_await sig;
}

cobalt::main co_main(int argc, char * argv[])
{

  signals::signal<void(int)> sig;
  auto p = await_signal(sig);
  sig(42);
  auto res = co_await p;
  assert(res == 42);
  co_return 0;
}