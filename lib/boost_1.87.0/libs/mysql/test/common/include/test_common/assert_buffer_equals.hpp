//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_ASSERT_BUFFER_EQUALS_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_ASSERT_BUFFER_EQUALS_HPP

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

namespace boost {
namespace mysql {
namespace test {

struct buffer_printer
{
    span<const std::uint8_t> buff;

    constexpr buffer_printer(span<const std::uint8_t> b) noexcept : buff(b) {}
};

std::ostream& operator<<(std::ostream& os, buffer_printer buff);
bool buffer_equals(span<const std::uint8_t> b1, span<const std::uint8_t> b2);

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
