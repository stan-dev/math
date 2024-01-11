//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/main.hpp>
#include <boost/outcome/coroutine_support.hpp>

using namespace boost;

outcome_v2::awaitables::lazy<int> lazy_func(int x)
{
  co_return x + 1;
}

outcome_v2::awaitables::eager<int> eager_func(int x)
{
  co_return x + 1;
}


cobalt::main co_main(int argc, char * argv[])
{
  [[maybe_unused]] auto lr = co_await lazy_func(10);

  assert(lr == 11);

  [[maybe_unused]] auto er = co_await eager_func(10);
  assert(er == 11);

  co_return 0;
}
