/* Boost.Flyweight test of concurrent_factory.
 *
 * Copyright 2024 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/flyweight for library home page.
 */

#include <boost/core/lightweight_test.hpp>
#include "test_concurrent_factory.hpp"

int main()
{
  test_concurrent_factory();
  return boost::report_errors();
}
