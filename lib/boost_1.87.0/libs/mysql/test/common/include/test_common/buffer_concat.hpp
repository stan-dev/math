//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_BUFFER_CONCAT_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_BUFFER_CONCAT_HPP

#include <boost/core/span.hpp>

#include <cstdint>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

inline std::vector<std::uint8_t> concat(std::vector<std::uint8_t> lhs, const std::vector<std::uint8_t>& rhs)
{
    lhs.insert(lhs.end(), rhs.begin(), rhs.end());
    return lhs;
}

class buffer_builder
{
    std::vector<std::uint8_t> buff_;

public:
    buffer_builder() = default;
    buffer_builder& add(span<const std::uint8_t> value)
    {
        buff_.insert(buff_.end(), value.begin(), value.end());
        return *this;
    }
    buffer_builder& add(const std::vector<std::uint8_t>& value)
    {
        return add(span<const std::uint8_t>(value));
    }
    std::vector<std::uint8_t> build() { return std::move(buff_); }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
