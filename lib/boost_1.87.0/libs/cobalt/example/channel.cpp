// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <boost/cobalt/channel.hpp>
#include <boost/cobalt/main.hpp>
#include <boost/cobalt/promise.hpp>

#include <iostream>

namespace cobalt = boost::cobalt;

// tag::channel_example[]
cobalt::promise<void> producer(cobalt::channel<int> & chan)
{
  for (int i = 0; i < 4; i++)
    co_await chan.write(i);

  chan.close();
}

cobalt::main co_main(int argc, char * argv[])
{
  cobalt::channel<int> c;

  auto p = producer(c);
  while (c.is_open())
    std::cout << co_await c.read() << std::endl;

  co_await p;
  co_return 0;
}
// end::channel_example[]
