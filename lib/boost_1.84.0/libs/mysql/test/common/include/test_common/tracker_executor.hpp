//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_TRACKER_EXECUTOR_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_TRACKER_EXECUTOR_HPP

#include <boost/asio/any_io_executor.hpp>

#include <cstddef>

namespace boost {
namespace mysql {
namespace test {

struct executor_info
{
    std::size_t num_posts{};
    std::size_t num_dispatches{};

    std::size_t total() const noexcept { return num_dispatches + num_posts; }
};

asio::any_io_executor create_tracker_executor(asio::any_io_executor inner, executor_info* tracked_values);
executor_info get_executor_info(const asio::any_io_executor& exec);

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
