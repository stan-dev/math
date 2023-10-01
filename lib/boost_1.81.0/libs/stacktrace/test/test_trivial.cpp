// Copyright Antony Polukhin, 2022-2022.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// See https://github.com/boostorg/stacktrace/issues/116

#include <boost/stacktrace/stacktrace.hpp>

#include <iostream>

int main() {
  std::cout << boost::stacktrace::stacktrace() << std::endl;
}
