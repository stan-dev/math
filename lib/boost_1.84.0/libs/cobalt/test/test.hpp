// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_COBALT_TEST2_HPP
#define BOOST_COBALT_TEST2_HPP

#include <boost/cobalt/task.hpp>
#include <boost/cobalt/run.hpp>
#include <boost/cobalt/spawn.hpp>

#include <boost/test/unit_test.hpp>


inline void test_run(boost::cobalt::task<void> (*func) ())
{
  using namespace boost;
#if !defined(BOOST_COBALT_NO_PMR)
  cobalt::pmr::unsynchronized_pool_resource res;
  cobalt::this_thread::set_default_resource(&res);
#endif
  {
    asio::io_context ctx;
    cobalt::this_thread::set_executor(ctx.get_executor());
    spawn(ctx, func(),
          +[](std::exception_ptr e)
          {
            BOOST_CHECK(e == nullptr);
          });
    std::size_t n;
    n = ctx.run();

    if (::getenv("BOOST_COBALT_BRUTE_FORCE"))
      while (n-- > 0)
      {
        ctx.restart();
        spawn(ctx, func(),
              +[](std::exception_ptr e)
              {
                BOOST_CHECK(e == nullptr);
              });
        for (std::size_t i = n; i > 0; i--)
          ctx.run_one();
      }
  }
#if !defined(BOOST_COBALT_NO_PMR)
  cobalt::this_thread::set_default_resource(cobalt::pmr::get_default_resource());
#endif
}

// tag::test_case_macro[]
#define CO_TEST_CASE(Function)                                                                                     \
static ::boost::cobalt::task<void> Function##_impl();                                                              \
BOOST_AUTO_TEST_CASE(Function)                                                                                     \
{                                                                                                                  \
    test_run(&Function##_impl);                                                                                           \
}                                                                                                                  \
static ::boost::cobalt::task<void> Function##_impl()
// end::test_case_macro[]

struct stop
{
  bool await_ready() {return false;}
  template<typename Promise>
  void await_suspend(std::coroutine_handle<Promise> h) { boost::cobalt::detail::self_destroy(h); }
  void await_resume() {}
};

struct immediate
{
  int state = 0;
  immediate() = default;
  immediate(const immediate & i);
  bool await_ready() {BOOST_CHECK(state++ == 0  ); return true;}
  void await_suspend(std::coroutine_handle<>) { BOOST_REQUIRE(false); }
  void await_resume() {BOOST_CHECK(state++ == 1);}

  ~immediate()
  {
    if (state != 0)
      BOOST_CHECK(state == 2);
  }
};

struct immediate_bool
{
  int state = 0;

  bool await_ready() {BOOST_CHECK(state++ == 0); return false;}
  bool await_suspend(std::coroutine_handle<>) { BOOST_CHECK(state++ == 1); return false; }
  void await_resume() {BOOST_CHECK(state++ == 2);}

  ~immediate_bool()
  {
    if (state != 0)
      BOOST_CHECK(state == 3);
  }
};

struct immediate_handle
{
  int state = 0;

  bool await_ready() {BOOST_CHECK(state++ == 0); return false;}
  std::coroutine_handle<> await_suspend(std::coroutine_handle<> h) { BOOST_CHECK(state++ == 1); return h; }
  void await_resume() {BOOST_CHECK(state++ == 2);}

  ~immediate_handle()
  {
    if (state != 0)
      BOOST_CHECK(state == 3);
  }
};


struct posted
{
  int state = 0;

  bool await_ready() {BOOST_CHECK(state++ == 0); return false;}
  void await_suspend(std::coroutine_handle<> h)
  {
    BOOST_CHECK(state++ == 1);
    boost::asio::post(boost::cobalt::this_thread::get_executor(), h);
  }
  void await_resume() {BOOST_CHECK(state++ == 2);}
  ~posted()
  {
    if (state != 0)
      BOOST_CHECK(state == 3);
  }
};

struct posted_bool
{
  int state = 0;

  bool await_ready() {BOOST_CHECK(state++ == 0); return false;}
  bool await_suspend(std::coroutine_handle<> h)
  {
    BOOST_CHECK(state++ == 1);
    boost::asio::post(boost::cobalt::this_thread::get_executor(), h);
    return true;
  }
  void await_resume() {BOOST_CHECK(state++ == 2);}
  ~posted_bool()
  {
    if (state != 0)
      BOOST_CHECK(state == 3);
  }
};

struct posted_handle
{
  int state = 0;

  bool await_ready() {BOOST_CHECK(state++ == 0); return false;}
  std::coroutine_handle<> await_suspend(std::coroutine_handle<> h)
  {
    BOOST_CHECK(state++ == 1);
    return boost::cobalt::detail::post_coroutine(
        boost::cobalt::this_thread::get_executor(), h
        );
  }
  void await_resume() {BOOST_CHECK(state++ == 2);}
  ~posted_handle()
  {
    if (state != 0)
      BOOST_CHECK(state == 3);
  }
};

#endif //BOOST_COBALT_TEST2_HPP
