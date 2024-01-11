// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/this_coro.hpp>
#include <boost/cobalt/detail/util.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"

using namespace boost;

struct sig_helper
{
  asio::cancellation_signal sig;
};

struct coro_feature_tester
      : sig_helper,
        cobalt::promise_cancellation_base<>,
        cobalt::promise_throw_if_cancelled_base,
        cobalt::enable_await_allocator<coro_feature_tester>
{
  using cobalt::promise_cancellation_base<>::await_transform;
  using cobalt::promise_throw_if_cancelled_base::await_transform;
  using cobalt::enable_await_allocator<coro_feature_tester>::await_transform;

  std::suspend_never initial_suspend() {return {};}
  std::suspend_never final_suspend() noexcept {return {};}

  coro_feature_tester(coro_feature_tester * & ref)
    : cobalt::promise_cancellation_base<>(sig.slot())
  {
    ref = this;
  }
  void return_void() {}
  void unhandled_exception() {throw;}
  void get_return_object() {}

#if !defined(BOOST_COBALT_NO_PMR)
  cobalt::pmr::unsynchronized_pool_resource res;
  using allocator_type = cobalt::pmr::polymorphic_allocator<void>;
  allocator_type get_allocator() {return alloc;}
  allocator_type alloc{&res};
#endif
};


namespace std
{
template<>
struct coroutine_traits<void, coro_feature_tester*> //< don't do THIS at
{
  using promise_type = coro_feature_tester;
};
}

#define SELF_TEST_CASE(Function)                                                                                       \
static void Function##_impl(coro_feature_tester * this_);                                                              \
BOOST_AUTO_TEST_CASE(Function)                                                                                         \
{                                                                                                                      \
    Function##_impl(nullptr);                                                                                          \
}                                                                                                                      \
static void Function##_impl(coro_feature_tester * this_)

BOOST_AUTO_TEST_SUITE(this_coro);

SELF_TEST_CASE(promise_cancellation_base)
{
  BOOST_CHECK(!this_->cancelled());
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::none);

  this_->sig.emit(asio::cancellation_type::terminal);

  BOOST_CHECK(this_->cancelled() == asio::cancellation_type::terminal);
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::terminal);

  co_await cobalt::this_coro::reset_cancellation_state();

  BOOST_CHECK(!this_->cancelled());
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::none);

  this_->sig.emit(asio::cancellation_type::terminal);

  BOOST_CHECK(this_->cancelled() == asio::cancellation_type::terminal);
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::terminal);

  co_await cobalt::this_coro::reset_cancellation_state(asio::enable_partial_cancellation());

  BOOST_CHECK(!this_->cancelled());
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::none);

  this_->sig.emit(asio::cancellation_type::total);

  BOOST_CHECK(!this_->cancelled());
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::none);
  this_->sig.emit(asio::cancellation_type::all);
  BOOST_CHECK(this_->cancelled() == (asio::cancellation_type::terminal | asio::cancellation_type::partial));


  co_await cobalt::this_coro::reset_cancellation_state(
      asio::enable_partial_cancellation(),
      asio::enable_terminal_cancellation());

  BOOST_CHECK(!this_->cancelled());
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::none);

  this_->sig.emit(asio::cancellation_type::total);

  BOOST_CHECK(!this_->cancelled());
  BOOST_CHECK(this_->cancellation_state().cancelled() == asio::cancellation_type::none);
  this_->sig.emit(asio::cancellation_type::all);
  BOOST_CHECK(this_->cancelled() == (asio::cancellation_type::terminal | asio::cancellation_type::partial));
}


SELF_TEST_CASE(promise_throw_if_cancelled_base)
{
  BOOST_CHECK(co_await asio::this_coro::throw_if_cancelled());
  co_await asio::this_coro::throw_if_cancelled(false);
  BOOST_CHECK(!co_await asio::this_coro::throw_if_cancelled());
}

#if !defined(BOOST_COBALT_NO_PMR)
SELF_TEST_CASE(enable_await_allocator)
{
  BOOST_CHECK(this_->get_allocator() == co_await cobalt::this_coro::allocator);
}
#endif


BOOST_AUTO_TEST_SUITE_END();