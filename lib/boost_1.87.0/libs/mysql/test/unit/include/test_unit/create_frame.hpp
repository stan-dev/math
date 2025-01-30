//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_FRAME_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_FRAME_HPP

#include <boost/core/span.hpp>

#include <cstdint>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

std::vector<std::uint8_t> create_frame(std::uint8_t seqnum, span<const std::uint8_t> body);

inline std::vector<std::uint8_t> create_frame(std::uint8_t seqnum, const std::vector<std::uint8_t>& body)
{
    return create_frame(seqnum, boost::span<const std::uint8_t>(body));
}

inline std::vector<std::uint8_t> create_empty_frame(std::uint8_t seqnum)
{
    return create_frame(seqnum, boost::span<const std::uint8_t>());
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
