// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/asio/any_completion_handler.hpp>
#include <boost/asio/any_io_executor.hpp>
#include <boost/asio/async_result.hpp>
#include <boost/asio/consign.hpp>
#include <boost/asio/error.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/this_coro.hpp>

#include <chrono>
#include <memory>

#include "boost/cobalt/main.hpp"
#include "boost/cobalt/op.hpp"

#include "test.hpp"

namespace cobalt = boost::cobalt;

BOOST_AUTO_TEST_SUITE(any_completion_token);

void cobalt_sleep_impl(
    boost::asio::any_completion_handler<void(boost::system::error_code)> handler,
    boost::asio::any_io_executor ex,
    std::chrono::nanoseconds duration
)
{
  auto timer = std::make_shared<boost::asio::steady_timer>(ex, duration);
  timer->async_wait(boost::asio::consign(std::move(handler), timer));
}

template <typename CompletionToken>
inline auto cobalt_sleep(
    boost::asio::any_io_executor ex,
    std::chrono::nanoseconds duration,
    CompletionToken&& token
)
{
  return boost::asio::async_initiate<CompletionToken, void(boost::system::error_code)>(
      cobalt_sleep_impl,
      token,
      std::move(ex),
      duration
  );
}


CO_TEST_CASE(sleep_any_cpl_token)
{
  co_await cobalt_sleep(co_await cobalt::this_coro::executor, std::chrono::milliseconds(1), cobalt::use_op);
}

BOOST_AUTO_TEST_SUITE_END();
