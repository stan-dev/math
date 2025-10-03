// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/detail/util.hpp>

#include <boost/cobalt/this_thread.hpp>
#include <boost/cobalt/unique_handle.hpp>

#include <boost/asio/bind_allocator.hpp>
#include <boost/asio/post.hpp>

namespace boost::cobalt::detail
{

#if BOOST_COBALT_NO_SELF_DELETE

void self_destroy(std::coroutine_handle<void> h, const cobalt::executor & exec) noexcept
{
#if defined(BOOST_COBALT_NO_PMR)
  asio::post(exec, [del=unique_handle<void>(h.address())]() mutable {});
#else
  asio::post(exec,
              asio::bind_allocator(
                 this_thread::get_allocator(),
                 [del=unique_handle<void>(h.address())]() mutable
                 {
                 }));
#endif

}
#endif

}
