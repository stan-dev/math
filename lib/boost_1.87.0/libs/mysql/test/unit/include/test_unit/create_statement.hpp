//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_STATEMENT_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_STATEMENT_HPP

#include <boost/mysql/statement.hpp>

#include <boost/mysql/detail/access.hpp>

#include <cstdint>

namespace boost {
namespace mysql {
namespace test {

class statement_builder
{
    std::uint32_t id_{};
    std::uint16_t num_params_{};

public:
    statement_builder() = default;
    statement_builder& id(std::uint32_t v) noexcept
    {
        id_ = v;
        return *this;
    }
    statement_builder& num_params(std::uint16_t v) noexcept
    {
        num_params_ = v;
        return *this;
    }
    statement build() { return detail::access::construct<statement>(id_, num_params_); }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
