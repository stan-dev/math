//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_FF_CHARSET_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_FF_CHARSET_HPP

#include <boost/mysql/character_set.hpp>

namespace boost {
namespace mysql {
namespace test {

// A hypothetical character set with rules that may confuse formatting algorithms.
// Some MySQL charsets (e.g. gbk) contain ASCII-compatible continuation characters, like this one.
inline std::size_t ff_charset_next_char(boost::span<const unsigned char> r)
{
    if (r[0] == 0xff)  // 0xff marks a two-byte character
        return r.size() > 1u ? 2u : 0u;
    return 1u;
};
constexpr character_set ff_charset{"ff_charset", ff_charset_next_char};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
