// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/detail/handler.hpp>

#include <boost/test/unit_test.hpp>
#include "test.hpp"
#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/post.hpp>

using namespace boost;

struct dummy_promise
{
    using executor_type = boost::asio::any_io_executor;
    executor_type get_executor() const;

};

static_assert(boost::asio::detail::has_executor_type<dummy_promise>::value);


void test(boost::cobalt::completion_handler<> ch)
{
  boost::asio::post(std::move(ch));
}

BOOST_AUTO_TEST_SUITE(handler);

struct immediate_aw
{
  bool await_ready() {return false;}

  std::optional<std::tuple<>> result;
  cobalt::detail::completed_immediately_t completed_immediately;

  template<typename T>
  bool await_suspend(std::coroutine_handle<T> h)
  {
    cobalt::completion_handler<> ch{h, result,
#if !defined(BOOST_COBALT_NO_PMR)
                                   cobalt::this_thread::get_default_resource(),
#endif
                                   &completed_immediately};

    auto exec = asio::get_associated_immediate_executor(ch, h.promise().get_executor());
    completed_immediately = cobalt::detail::completed_immediately_t::initiating;
    asio::dispatch(exec, std::move(ch));

    BOOST_CHECK(result);
    BOOST_CHECK(completed_immediately == cobalt::detail::completed_immediately_t::yes);

    return completed_immediately != cobalt::detail::completed_immediately_t::yes;
  }

  void await_resume()
  {
    BOOST_CHECK(completed_immediately != cobalt::detail::completed_immediately_t::no);
    BOOST_CHECK(result);
  }
};

#if !defined(BOOST_COBALT_USE_IO_CONTEXT)

struct non_immediate_aw
{
  bool await_ready() {return false;}

  std::optional<std::tuple<>> result;
  cobalt::detail::completed_immediately_t completed_immediately;
  cobalt::detail::sbo_resource res;

  template<typename T>
  bool await_suspend(std::coroutine_handle<T> h)
  {
    cobalt::completion_handler<> ch{h, result, &res, &completed_immediately};

    auto exec = asio::get_associated_immediate_executor(ch, h.promise().get_executor());
    asio::dispatch(exec,
                   asio::deferred(
                       [exec = h.promise().get_executor()]()
                       {
                         return asio::post(exec, asio::deferred);
                       }))(std::move(ch));

    BOOST_CHECK(!result);
    BOOST_CHECK(completed_immediately != cobalt::detail::completed_immediately_t::yes);

    return completed_immediately != cobalt::detail::completed_immediately_t::yes;
  }

  void await_resume()
  {
    BOOST_CHECK(completed_immediately != cobalt::detail::completed_immediately_t::yes);
    BOOST_CHECK(result);
  }
};




CO_TEST_CASE(immediate_completion)
{
  co_await immediate_aw{};
  co_await non_immediate_aw{};
}

#endif

BOOST_AUTO_TEST_CASE(immediate_executor_initiating)
{
  asio::io_context ctx;
  cobalt::detail::completed_immediately_t completed_immediately = cobalt::detail::completed_immediately_t::initiating;
  cobalt::detail::completion_handler_noop_executor chh{ctx.get_executor(), &completed_immediately};
  bool called = false;

  asio::dispatch(chh, [&] { called = true; });
  BOOST_CHECK(called);
  BOOST_CHECK(completed_immediately == cobalt::detail::completed_immediately_t::initiating);
}

BOOST_AUTO_TEST_CASE(immediate_executor_maybe)
{
  asio::io_context ctx;
  cobalt::detail::completed_immediately_t completed_immediately = cobalt::detail::completed_immediately_t::initiating;
  cobalt::detail::completion_handler_noop_executor chh{ctx.get_executor(), &completed_immediately};
  bool called = false;

  completed_immediately = cobalt::detail::completed_immediately_t::maybe;
  asio::dispatch(chh, [&] { called = true; completed_immediately = cobalt::detail::completed_immediately_t::yes; });
  BOOST_CHECK(called);
  BOOST_CHECK(completed_immediately == cobalt::detail::completed_immediately_t::yes);
}


BOOST_AUTO_TEST_CASE(immediate_executor_maybe_not)
{
  asio::io_context ctx;
  cobalt::detail::completed_immediately_t completed_immediately = cobalt::detail::completed_immediately_t::initiating;
  cobalt::detail::completion_handler_noop_executor chh{ctx.get_executor(), &completed_immediately};
  bool called = false;

  completed_immediately = cobalt::detail::completed_immediately_t::maybe;
  asio::dispatch(chh, [&] { called = true; });
  BOOST_CHECK(called);
  BOOST_CHECK(completed_immediately == cobalt::detail::completed_immediately_t::initiating);
}

BOOST_AUTO_TEST_CASE(immediate_executor_no)
{
  asio::io_context ctx;
  cobalt::detail::completed_immediately_t completed_immediately = cobalt::detail::completed_immediately_t::initiating;
  cobalt::detail::completion_handler_noop_executor chh{ctx.get_executor(), &completed_immediately};
  bool called = false;

  completed_immediately = cobalt::detail::completed_immediately_t::no;
  asio::dispatch(chh, [&] { called = true; });
  BOOST_CHECK(!called);
  BOOST_CHECK(completed_immediately == cobalt::detail::completed_immediately_t::no);
  BOOST_CHECK(ctx.run() == 1u);
  BOOST_CHECK(called);
}

BOOST_AUTO_TEST_SUITE_END();
