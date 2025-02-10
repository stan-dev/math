// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <boost/cobalt.hpp>
#include <boost/asio.hpp>

#if defined(BOOST_COBALT_BENCH_WITH_CONTEXT)
#include <boost/asio/spawn.hpp>
#endif

using namespace boost;
constexpr std::size_t n = 50'000'000ull;

cobalt::task<void> atest()
{
  for (std::size_t i = 0u; i < n; i++)
    co_await asio::post(cobalt::use_op);
}

asio::awaitable<void> awtest()
{
  for (std::size_t i = 0u; i < n; i++)
    co_await asio::post(asio::deferred);
}

#if defined(BOOST_COBALT_BENCH_WITH_CONTEXT)

void stest(asio::yield_context ctx)
{
  for (std::size_t i = 0u; i < n; i++)
    asio::post(ctx);
}

#endif


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

#if defined(BOOST_COBALT_BENCH_WITH_CONTEXT)
  {
    auto start = std::chrono::steady_clock::now();
    asio::io_context ctx{BOOST_ASIO_CONCURRENCY_HINT_1};
    asio::spawn(ctx, stest, asio::detached);
    ctx.run();
    auto end = std::chrono::steady_clock::now();
    printf("stackful : %ld ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  }
#endif

  return 0;
}