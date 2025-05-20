// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/error.hpp>

namespace boost::cobalt
{
system::error_category & cobalt_category()
{
  static cobalt_category_t cat;
  return cat;
}

system::error_code make_error_code(error e)
{
  return system::error_code(static_cast<int>(e), cobalt_category());
}

}

