//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_ASIO_NO_DEPRECATED 1
#include <boost/process.hpp>
int main() {}

#if defined(BOOST_POSIX_API)
#include <boost/process/posix.hpp>
#else
#include <boost/process/windows.hpp>
#endif