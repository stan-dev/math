//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/this_thread.hpp>
#include <boost/asio/any_io_executor.hpp>

#include <boost/asio/executor.hpp>



#include <optional>

namespace boost::cobalt::this_thread
{

namespace detail
{
#if !defined(BOOST_COBALT_NO_PMR)
thread_local pmr::memory_resource *default_coro_memory_resource = pmr::get_default_resource();
#endif

thread_local std::optional<executor> executor;
}

#if !defined(BOOST_COBALT_NO_PMR)
pmr::memory_resource* get_default_resource() noexcept
{
  return detail::default_coro_memory_resource;
}

pmr::memory_resource* set_default_resource(pmr::memory_resource* r) noexcept
{
  auto pre = get_default_resource();
  detail::default_coro_memory_resource = r;
  return pre;
}

pmr::polymorphic_allocator<void> get_allocator()
{
  return pmr::polymorphic_allocator<void>{get_default_resource()};
}
#endif

bool has_executor()
{
  return detail::executor.has_value();
}

executor & get_executor(const boost::source_location & loc)
{
  if (!detail::executor)
    throw_exception(asio::bad_executor(), loc);
  return *detail::executor;
}

struct this_thread_service : asio::detail::execution_context_service_base<this_thread_service>
{
  this_thread_service(asio::execution_context & ctx)
      : asio::detail::execution_context_service_base<this_thread_service>(ctx)
  {
  }


  void shutdown() override
  {
    if (detail::executor && (&detail::executor->context() == &this->context()))
      detail::executor.reset();
  }
};

void set_executor(executor exec) noexcept
{
  detail::executor = std::move(exec);
  asio::use_service<this_thread_service>(detail::executor->context());
}
}

namespace boost::cobalt::detail
{

#if defined(BOOST_COBALT_CUSTOM_EXECUTOR) || defined(BOOST_COBALT_USE_IO_CONTEXT)
executor
extract_executor(asio::any_io_executor exec)
{
  auto t = exec.target<executor>();
  if (t == nullptr)
    throw_exception(asio::bad_executor());

  return *t;
}
#endif
}
