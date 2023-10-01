//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_BUFFER_CONCAT_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_BUFFER_CONCAT_HPP

#include <boost/core/span.hpp>

#include <cstdint>
#include <cstring>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

inline void concat(std::vector<std::uint8_t>& lhs, const void* buff, std::size_t size)
{
    if (size)
    {
        auto current_size = lhs.size();
        lhs.resize(current_size + size);
        std::memcpy(lhs.data() + current_size, buff, size);
    }
}

inline void concat(std::vector<std::uint8_t>& lhs, const std::vector<uint8_t>& rhs)
{
    concat(lhs, rhs.data(), rhs.size());
}

inline std::vector<std::uint8_t> concat_copy(
    std::vector<std::uint8_t> lhs,
    const std::vector<std::uint8_t>& rhs
)
{
    concat(lhs, rhs);
    return lhs;
}

class buffer_builder
{
    std::vector<std::uint8_t> buff_;

public:
    buffer_builder() = default;
    buffer_builder& add(span<const std::uint8_t> value)
    {
        concat(buff_, value.data(), value.size());
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
