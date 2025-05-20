//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_CHECK_META_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_CHECK_META_HPP

// This is a lighter check than integ tests' metadata_validator

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/test/unit_test.hpp>

#include <initializer_list>

namespace boost {
namespace mysql {
namespace test {

inline void check_meta(metadata_collection_view meta, const std::vector<column_type>& expected_types)
{
    BOOST_TEST_REQUIRE(meta.size() == expected_types.size());
    for (std::size_t i = 0; i < meta.size(); ++i)
    {
        BOOST_TEST(meta[i].type() == expected_types[i]);
    }
}

inline void check_meta(
    metadata_collection_view meta,
    const std::vector<std::pair<column_type, string_view>>& expected
)
{
    BOOST_TEST_REQUIRE(meta.size() == expected.size());
    for (std::size_t i = 0; i < meta.size(); ++i)
    {
        BOOST_TEST(meta[i].type() == expected[i].first);
        BOOST_TEST(meta[i].type() == expected[i].first);
    }
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
