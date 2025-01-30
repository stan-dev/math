//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_COLDEF_FRAME_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_COLDEF_FRAME_HPP

#include <boost/mysql/detail/coldef_view.hpp>

#include "test_unit/create_frame.hpp"

namespace boost {
namespace mysql {
namespace test {

std::vector<std::uint8_t> create_coldef_body(const detail::coldef_view& pack);

inline std::vector<std::uint8_t> create_coldef_frame(std::uint8_t seqnum, const detail::coldef_view& coldef)
{
    return create_frame(seqnum, create_coldef_body(coldef));
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
