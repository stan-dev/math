//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_VALIDATE_STRING_CONTAINS_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_VALIDATE_STRING_CONTAINS_HPP

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

inline void validate_string_contains(std::string value, const std::vector<std::string>& to_check)
{
    std::transform(value.begin(), value.end(), value.begin(), [](char c) {
        return static_cast<char>(tolower(c));
    });
    for (const auto& elm : to_check)
    {
        BOOST_TEST(
            value.find(elm) != std::string::npos,
            "Substring '" << elm << "' not found in '" << value << "'"
        );
    }
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
