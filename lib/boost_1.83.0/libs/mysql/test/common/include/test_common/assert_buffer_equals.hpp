//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_ASSERT_BUFFER_EQUALS_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_ASSERT_BUFFER_EQUALS_HPP

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <cstring>
#include <iomanip>

namespace boost {
namespace mysql {
namespace test {

struct buffer_printer
{
    span<const std::uint8_t> buff;

    constexpr buffer_printer(span<const std::uint8_t> b) noexcept : buff(b) {}
};

inline std::ostream& operator<<(std::ostream& os, buffer_printer buff)
{
    os << std::setfill('0') << std::hex << "{ ";
    for (std::size_t i = 0; i < buff.buff.size(); ++i)
    {
        os << "0x" << std::setw(2) << static_cast<int>(buff.buff.data()[i]) << ", ";
    }
    return os << "}";
}

inline bool buffer_equals(span<const std::uint8_t> b1, span<const std::uint8_t> b2) noexcept
{
    // If any of the buffers are empty (data() == nullptr), prevent
    // calling memcmp (UB)
    if (b1.size() == 0 || b2.size() == 0)
        return b1.size() == 0 && b2.size() == 0;

    if (b1.size() != b2.size())
        return false;

    return ::std::memcmp(b1.data(), b2.data(), b1.size()) == 0;
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#define BOOST_MYSQL_ASSERT_BUFFER_EQUALS(b1, b2)                                            \
    BOOST_TEST(                                                                             \
        ::boost::mysql::test::buffer_equals(b1, b2),                                        \
        #b1 " != " #b2 ": \nlhs: " << ::boost::mysql::test::buffer_printer(b1)              \
                                   << "\nrhs: " << ::boost::mysql::test::buffer_printer(b2) \
    )

#endif
