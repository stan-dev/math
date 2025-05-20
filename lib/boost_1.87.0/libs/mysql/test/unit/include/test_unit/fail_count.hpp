//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_FAIL_COUNT_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_FAIL_COUNT_HPP

#include <boost/mysql/error_code.hpp>

namespace boost {
namespace mysql {
namespace test {

// Inspired by Beast's fail count
class fail_count
{
    std::size_t fail_after_;
    std::size_t num_calls_{0};
    error_code err_;

public:
    static constexpr std::size_t never_fail = std::size_t(-1);
    explicit fail_count(
        std::size_t fail_after = never_fail,
        error_code err = make_error_code(std::errc::io_error)
    ) noexcept
        : fail_after_(fail_after), err_(err)
    {
    }
    error_code maybe_fail() noexcept { return num_calls_++ >= fail_after_ ? err_ : error_code(); }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
