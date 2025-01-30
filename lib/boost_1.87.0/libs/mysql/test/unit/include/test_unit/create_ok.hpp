//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_OK_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_OK_HPP

#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/flags.hpp>
#include <boost/mysql/detail/ok_view.hpp>

namespace boost {
namespace mysql {
namespace test {

class ok_builder
{
    detail::ok_view ok_{};

    void flag(std::uint16_t f, bool value) noexcept
    {
        if (value)
            ok_.status_flags |= f;
        else
            ok_.status_flags &= ~f;
    }

public:
    ok_builder() = default;
    ok_builder& affected_rows(std::uint64_t v) noexcept
    {
        ok_.affected_rows = v;
        return *this;
    }
    ok_builder& last_insert_id(std::uint64_t v) noexcept
    {
        ok_.last_insert_id = v;
        return *this;
    }
    ok_builder& warnings(std::uint16_t v) noexcept
    {
        ok_.warnings = v;
        return *this;
    }
    ok_builder& flags(std::uint16_t v) noexcept
    {
        ok_.status_flags = v;
        return *this;
    }
    ok_builder& more_results(bool v) noexcept
    {
        flag(detail::status_flags::more_results, v);
        return *this;
    }
    ok_builder& out_params(bool v) noexcept
    {
        flag(detail::status_flags::out_params, v);
        return *this;
    }
    ok_builder& no_backslash_escapes(bool v) noexcept
    {
        flag(detail::status_flags::no_backslash_escapes, v);
        return *this;
    }
    ok_builder& info(string_view v) noexcept
    {
        ok_.info = v;
        return *this;
    }
    detail::ok_view build() const noexcept { return ok_; }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
