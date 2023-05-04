//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_BUFFER_CONCAT_HPP
#define BOOST_MYSQL_TEST_COMMON_BUFFER_CONCAT_HPP

#include <cstdint>
#include <cstring>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

inline void concat(std::vector<std::uint8_t>& lhs, const void* buff, std::size_t size)
{
    auto current_size = lhs.size();
    lhs.resize(current_size + size);
    std::memcpy(lhs.data() + current_size, buff, size);
}

inline void concat(std::vector<std::uint8_t>& lhs, const std::vector<uint8_t>& rhs)
{
    concat(lhs, rhs.data(), rhs.size());
}

inline std::vector<std::uint8_t> concat_copy(
    std::vector<uint8_t> lhs,
    const std::vector<uint8_t>& rhs
)
{
    concat(lhs, rhs);
    return lhs;
}

inline std::vector<std::uint8_t> concat_copy(
    std::vector<uint8_t> lhs,
    const std::vector<uint8_t>& rhs,
    const std::vector<uint8_t>& rhs2
)
{
    concat(lhs, rhs);
    concat(lhs, rhs2);
    return lhs;
}

inline std::vector<std::uint8_t> concat_copy(
    std::vector<uint8_t> lhs,
    const std::vector<uint8_t>& rhs,
    const std::vector<uint8_t>& rhs2,
    const std::vector<uint8_t>& rhs3
)
{
    concat(lhs, rhs);
    concat(lhs, rhs2);
    concat(lhs, rhs3);
    return lhs;
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
