//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_MOCK_MESSAGE_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_MOCK_MESSAGE_HPP

#include <boost/mysql/impl/internal/protocol/impl/serialization_context.hpp>

#include <boost/core/span.hpp>

#include <cstdint>

namespace boost {
namespace mysql {
namespace test {

struct mock_message
{
    span<const std::uint8_t> data;

    void serialize(detail::serialization_context& ctx) const { ctx.add(data); }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
