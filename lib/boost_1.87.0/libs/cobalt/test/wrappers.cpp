// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/detail/wrapper.hpp>

#include <boost/asio/detached.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/bind_allocator.hpp>


#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(wrappers);

BOOST_AUTO_TEST_CASE(regular)
{
    boost::asio::io_context ctx;
    boost::cobalt::this_thread::set_executor(ctx.get_executor());
    bool ran = false;

#if !defined(BOOST_COBALT_NO_PMR)
    char buf[512];
    boost::cobalt::pmr::monotonic_buffer_resource res{buf, 512};
    struct completion
    {
      bool & ran;
      using allocator_type = boost::cobalt::pmr::polymorphic_allocator<void>;
      allocator_type get_allocator() const { return alloc; }
      boost::cobalt::pmr::polymorphic_allocator<void> alloc;
      void operator()()
      {
        ran = true;
      }
    };

    auto p = boost::cobalt::detail::post_coroutine(ctx.get_executor(),
                                              completion{ran, boost::cobalt::pmr::polymorphic_allocator<void>(&res)}
                                          );
#else
    auto p = boost::cobalt::detail::post_coroutine(ctx.get_executor(), [&]{ran = true;});
#endif
    BOOST_CHECK(p);
    BOOST_CHECK(!ran);
    p.resume();
    BOOST_CHECK(!ran);
    ctx.run();
    BOOST_CHECK(ran);
}

BOOST_AUTO_TEST_CASE(expire)
{

  boost::asio::io_context ct2;
  boost::cobalt::this_thread::set_executor(ct2.get_executor());
  auto h = boost::cobalt::detail::post_coroutine(ct2.get_executor(), boost::asio::detached);
  boost::cobalt::detail::self_destroy(h);
}


BOOST_AUTO_TEST_SUITE_END();