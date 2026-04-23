// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <boost/cobalt.hpp>
#include <boost/asio.hpp>
#include <boost/asio/experimental/channel.hpp>
#include <boost/asio/experimental/parallel_group.hpp>
#include <boost/asio/experimental/awaitable_operators.hpp>

#if defined(BOOST_COBALT_BENCH_WITH_CONTEXT)
#include <boost/asio/spawn.hpp>
#endif

using namespace boost;
constexpr std::size_t n = 3'000'000ull;

cobalt::task<void> atest()
{
  asio::experimental::channel<void(system::error_code)> chan{co_await cobalt::this_coro::executor, 0u};
  for (std::size_t i = 0u; i < n; i++)
    co_await cobalt::gather(
              chan.async_send(system::error_code{}, cobalt::use_task),
              chan.async_receive(cobalt::use_task));

}

asio::awaitable<void> awtest()
{
  asio::experimental::channel<void(system::error_code)> chan{co_await cobalt::this_coro::executor, 0u};
  using boost::asio::experimental::awaitable_operators::operator&&;
  for (std::size_t i = 0u; i < n; i++)
    co_await (
        chan.async_send(system::error_code{}, asio::use_awaitable)
        &&
        chan.async_receive(asio::use_awaitable));
}

int main(int argc, char * argv[])
{
  {
    auto start = std::chrono::steady_clock::now();
    cobalt::run(atest());
    auto end = std::chrono::steady_clock::now();
    printf("cobalt    : %ld ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  }

  {
    auto start = std::chrono::steady_clock::now();
    asio::io_context ctx{BOOST_ASIO_CONCURRENCY_HINT_1};
    asio::co_spawn(ctx, awtest(), asio::detached);
    ctx.run();
    auto end = std::chrono::steady_clock::now();
    printf("awaitable: %ld ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  }

  return 0;
}