//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_META_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_META_HPP

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/flags.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

class meta_builder
{
    detail::coldef_view coldef_{};

public:
    meta_builder()
    {
        coldef_.collation_id = 33;  // utf8_general_ci
        coldef_.type = column_type::enum_;
    }
    meta_builder& database(string_view v) noexcept
    {
        coldef_.database = v;
        return *this;
    }
    meta_builder& table(string_view v) noexcept
    {
        coldef_.table = v;
        return *this;
    }
    meta_builder& org_table(string_view v) noexcept
    {
        coldef_.org_table = v;
        return *this;
    }
    meta_builder& name(string_view v) noexcept
    {
        coldef_.name = v;
        return *this;
    }
    meta_builder& org_name(string_view v) noexcept
    {
        coldef_.org_name = v;
        return *this;
    }
    meta_builder& collation_id(std::uint16_t v) noexcept
    {
        coldef_.collation_id = v;
        return *this;
    }
    meta_builder& column_length(std::uint32_t v) noexcept
    {
        coldef_.column_length = v;
        return *this;
    }
    meta_builder& type(column_type v) noexcept
    {
        coldef_.type = v;
        return *this;
    }
    meta_builder& flags(std::uint16_t v) noexcept
    {
        coldef_.flags = v;
        return *this;
    }
    meta_builder& unsigned_flag(bool v) noexcept
    {
        if (v)
            coldef_.flags |= detail::column_flags::unsigned_;
        else
            coldef_.flags &= ~detail::column_flags::unsigned_;
        return *this;
    }
    meta_builder& nullable(bool v) noexcept
    {
        if (v)
            coldef_.flags &= ~detail::column_flags::not_null;
        else
            coldef_.flags |= detail::column_flags::not_null;
        return *this;
    }
    meta_builder& zerofill(bool v) noexcept
    {
        if (v)
            coldef_.flags &= ~detail::column_flags::zerofill;
        else
            coldef_.flags |= detail::column_flags::zerofill;
        return *this;
    }
    meta_builder& decimals(std::uint8_t v) noexcept
    {
        coldef_.decimals = v;
        return *this;
    }
    metadata build() const { return detail::access::construct<metadata>(coldef_, true); }
    detail::coldef_view build_coldef() const noexcept { return coldef_; }
};

inline metadata create_meta(column_type type) { return meta_builder().type(type).build(); }
inline std::vector<metadata> create_metas(const std::vector<column_type>& types)
{
    std::vector<metadata> res;
    res.reserve(types.size());
    for (column_type t : types)
        res.push_back(create_meta(t));
    return res;
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
