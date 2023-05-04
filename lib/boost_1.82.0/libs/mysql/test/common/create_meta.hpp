//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_CREATE_META_HPP
#define BOOST_MYSQL_TEST_COMMON_CREATE_META_HPP

#include <boost/mysql/metadata.hpp>
#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/protocol/common_messages.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>

#include <cstdint>

namespace boost {
namespace mysql {
namespace test {

inline metadata create_meta(const detail::column_definition_packet& coldef, bool copy_strings)
{
    return detail::metadata_access::construct(coldef, copy_strings);
}

inline metadata create_meta(
    detail::protocol_field_type type,
    std::uint16_t flags = 0,
    std::uint8_t decimals = 0,
    std::uint16_t collation = mysql_collations::utf8mb4_general_ci
)
{
    detail::column_definition_packet coldef{};
    coldef.type = type;
    coldef.flags = flags;
    coldef.decimals = decimals;
    coldef.character_set = collation;
    return create_meta(coldef, true);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
