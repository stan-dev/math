//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_CHANNEL_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_CHANNEL_HPP

#include <boost/mysql/detail/any_stream_impl.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>

#include <cstddef>
#include <memory>

#include "test_unit/test_stream.hpp"

namespace boost {
namespace mysql {
namespace test {

inline detail::channel create_channel(std::size_t buffer_size = 1024)
{
    return detail::channel(buffer_size, std::unique_ptr<test_any_stream>(new test_any_stream));
}

inline test_stream& get_stream(detail::channel& chan) noexcept
{
    return detail::cast<test_stream>(chan.stream());
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
