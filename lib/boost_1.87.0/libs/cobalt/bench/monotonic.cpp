// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <boost/cobalt.hpp>
#include <boost/cobalt/detail/monotonic_resource.hpp>
#include <boost/cobalt/detail/sbo_resource.hpp>
#include <boost/asio.hpp>
#include <boost/asio/yield.hpp>


#include <boost/core/demangle.hpp>

using namespace boost;
constexpr std::size_t n = 2'000'000ull;


/* The SBO optimization is for structures allocations,
 * i.e. in a LIFO structure.
*/
struct my_composed_op : asio::coroutine
{
  std::shared_ptr<char[]> res1, res2;

  template<typename Self>
  void operator()(Self && self)
  {
    reenter(this)
    {
      yield asio::post(std::move(self));
      res1 = std::allocate_shared<char[]>(self.get_allocator(), 512);
      yield asio::post(std::move(self));
      res2 = std::allocate_shared<char[]>(self.get_allocator(), 256);
      yield asio::post(std::move(self));
      res2.reset(); // tests reusing memory
      res2 = std::allocate_shared<char[]>(self.get_allocator(), 256);
      yield asio::post(std::move(self));
      res2.reset(); // tests reusing memory
      res2 = std::allocate_shared<char[]>(self.get_allocator(), 256);
      yield asio::post(std::move(self));
      res2.reset(); // tests reusing memory
      res2 = std::allocate_shared<char[]>(self.get_allocator(), 256);
      self.complete();
    }
  }
};

template<typename Token>
auto my_composed(asio::io_context & ctx, Token && token)
{
  return asio::cobalt_compose<Token, void()>(my_composed_op{}, token, ctx.get_executor());
}


struct std_test
{
  asio::io_context & ctx;
  std::size_t i = 0u;
  void operator()()
  {
    if (i ++ < n)
      my_composed(ctx, std::move(*this));
  }
};

alignas(std::max_align_t) char buf[1024];
cobalt::detail::monotonic_resource res{buf, sizeof(buf)};

struct mono_test
{
  asio::io_context & ctx;
  std::size_t i = 0u;

  using allocator_type = cobalt::detail::monotonic_allocator<void>;
  allocator_type get_allocator() const { return cobalt::detail::monotonic_allocator<void>{&res}; }

  void operator()()
  {
    res.release();  // this is a thing that wouldn't happen in a comopsed op
    if (i ++ < n)
      my_composed(ctx, std::move(*this));
  }
};

cobalt::pmr::monotonic_buffer_resource pmr_res{buf, sizeof(buf)};

struct pmr_test
{
  asio::io_context & ctx;
  std::size_t i = 0u;

  using allocator_type = cobalt::pmr::polymorphic_allocator<void>;
  allocator_type get_allocator() const { return cobalt::pmr::polymorphic_allocator<void>{&pmr_res}; }

  void operator()()
  {
    pmr_res.release(); // this is a thing that wouldn't happen in a composed op
    if (i ++ < n)
      my_composed(ctx, std::move(*this));
  }
};


cobalt::detail::sbo_resource sbo_res{buf, sizeof(buf)};

struct sbo_test
{
  asio::io_context & ctx;
  std::size_t i = 0u;

  using allocator_type = cobalt::detail::sbo_allocator<void>;
  allocator_type get_allocator() const { return cobalt::detail::sbo_allocator<void>{&sbo_res}; }

  void operator()()
  {
    if (i ++ < n)
      my_composed(ctx, std::move(*this));
  }
};




int main(int argc, char * argv[])
{

  {
    auto start = std::chrono::steady_clock::now();
    asio::io_context ctx{BOOST_ASIO_CONCURRENCY_HINT_1};
    std_test{ctx}();
    ctx.run();
    auto end = std::chrono::steady_clock::now();
    printf("std::allocator  : %ld ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  }

  {
    auto start = std::chrono::steady_clock::now();
    asio::io_context ctx{BOOST_ASIO_CONCURRENCY_HINT_1};
    mono_test{ctx}();
    ctx.run();
    auto end = std::chrono::steady_clock::now();
    printf("cobalt::monotonic: %ld ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  }

  {
    auto start = std::chrono::steady_clock::now();
    asio::io_context ctx{BOOST_ASIO_CONCURRENCY_HINT_1};
    pmr_test{ctx}();
    ctx.run();
    auto end = std::chrono::steady_clock::now();
    printf("pmr::monotonic: %ld ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  }


  {
    auto start = std::chrono::steady_clock::now();
    asio::io_context ctx{BOOST_ASIO_CONCURRENCY_HINT_1};
    sbo_test{ctx}();
    ctx.run();
    auto end = std::chrono::steady_clock::now();
    printf("cobalt::sbo: %ld ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
  }


  return 0;
}