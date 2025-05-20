//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/detail/exception.hpp>
#include <boost/cobalt/error.hpp>

namespace boost::cobalt::detail
{

std::exception_ptr moved_from_exception()
{
  static auto ep = std::make_exception_ptr(system::error_code(
      error::moved_from
      ));
  return ep;
}

std::exception_ptr detached_exception()
{

  static auto ep = std::make_exception_ptr(system::error_code(
        error::detached
      ));
  return ep;
}

std::exception_ptr completed_unexpected()
{

  static auto ep = std::make_exception_ptr(system::error_code(
          error::completed_unexpected
      ));
  return ep;
}

std::exception_ptr wait_not_ready()
{
  static auto ep = std::make_exception_ptr(system::error_code(
          error::wait_not_ready
      ));
  return ep;
}

std::exception_ptr already_awaited()
{
  static auto ep = std::make_exception_ptr(system::error_code(
        error::already_awaited
    ));
  return ep;
}


std::exception_ptr allocation_failed()
{
  static auto ep = std::make_exception_ptr(system::error_code(
      error::already_awaited
  ));
  return ep;
}


}
