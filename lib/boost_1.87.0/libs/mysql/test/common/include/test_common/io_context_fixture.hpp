//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_IO_CONTEXT_FIXTURE_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_IO_CONTEXT_FIXTURE_HPP

#include <boost/asio/io_context.hpp>

namespace boost {
namespace mysql {
namespace test {

struct io_context_fixture
{
    asio::io_context ctx;

    // Checks that we effectively run out of work
    ~io_context_fixture();
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
