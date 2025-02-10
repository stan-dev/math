//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/metadata.hpp>

#include <boost/mysql/detail/typing/pos_map.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>

#include "test_common/create_basic.hpp"
#include "test_unit/create_meta.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::span;
using boost::mysql::detail::map_field_view;
using boost::mysql::detail::map_metadata;
using boost::mysql::detail::name_table_t;
using boost::mysql::detail::pos_absent;
using boost::mysql::detail::pos_map_add_field;
using boost::mysql::detail::pos_map_reset;

BOOST_AUTO_TEST_SUITE(test_post_map)

BOOST_AUTO_TEST_CASE(reset_empty)
{
    span<std::size_t> map{};
    BOOST_CHECK_NO_THROW(pos_map_reset(map));
}

BOOST_AUTO_TEST_CASE(reset_nonempty)
{
    std::array<std::size_t, 4> storage{
        {42, 43, 44, 45}
    };
    span<std::size_t> map(storage.data(), 3);

    pos_map_reset(map);
    BOOST_TEST(map[0] == pos_absent);
    BOOST_TEST(map[1] == pos_absent);
    BOOST_TEST(map[2] == pos_absent);
    BOOST_TEST(map[3] == 45u);  // didn't modify any extra storage
}

BOOST_AUTO_TEST_CASE(add_field_empty)
{
    span<std::size_t> map{};
    name_table_t name_table{};
    BOOST_CHECK_NO_THROW(pos_map_add_field(map, name_table, 0, "f1"));
}

BOOST_AUTO_TEST_CASE(add_field_unnamed)
{
    // Setup
    std::array<std::size_t, 4> map{
        {42, 43, 44}
    };
    name_table_t name_table{};
    pos_map_reset(map);

    // Add first field
    pos_map_add_field(map, name_table, 0, "f1");
    BOOST_TEST(map[0] == 0u);
    BOOST_TEST(map[1] == pos_absent);
    BOOST_TEST(map[2] == pos_absent);

    // Add second field
    pos_map_add_field(map, name_table, 1, "f2");
    BOOST_TEST(map[0] == 0u);
    BOOST_TEST(map[1] == 1u);
    BOOST_TEST(map[2] == pos_absent);

    // Add third field
    pos_map_add_field(map, name_table, 2, "f3");
    BOOST_TEST(map[0] == 0u);
    BOOST_TEST(map[1] == 1u);
    BOOST_TEST(map[2] == 2u);

    // Any further trailing fields are discarded
    BOOST_CHECK_NO_THROW(pos_map_add_field(map, name_table, 3, "f3"));
    BOOST_CHECK_NO_THROW(pos_map_add_field(map, name_table, 4, "f4"));
    BOOST_TEST(map[0] == 0u);
    BOOST_TEST(map[1] == 1u);
    BOOST_TEST(map[2] == 2u);
}

BOOST_AUTO_TEST_CASE(add_field_named)
{
    // Setup
    const string_view name_table[] = {"f1", "f2", "f3", "f4"};
    std::array<std::size_t, 4> map{{}};
    pos_map_reset(map);

    // Add first field
    pos_map_add_field(map, name_table, 0, "f2");
    BOOST_TEST(map[0] == pos_absent);
    BOOST_TEST(map[1] == 0u);
    BOOST_TEST(map[2] == pos_absent);
    BOOST_TEST(map[3] == pos_absent);

    // Add second field
    pos_map_add_field(map, name_table, 1, "f4");
    BOOST_TEST(map[0] == pos_absent);
    BOOST_TEST(map[1] == 0u);
    BOOST_TEST(map[2] == pos_absent);
    BOOST_TEST(map[3] == 1u);

    // Add a non-existing field
    pos_map_add_field(map, name_table, 2, "fnonexistent");
    BOOST_TEST(map[0] == pos_absent);
    BOOST_TEST(map[1] == 0u);
    BOOST_TEST(map[2] == pos_absent);
    BOOST_TEST(map[3] == 1u);

    // Add third field
    pos_map_add_field(map, name_table, 3, "f1");
    BOOST_TEST(map[0] == 3u);
    BOOST_TEST(map[1] == 0u);
    BOOST_TEST(map[2] == pos_absent);
    BOOST_TEST(map[3] == 1u);
}

BOOST_AUTO_TEST_CASE(map_metadata_)
{
    const std::array<std::size_t, 3> map{
        {1, 0, 2}
    };
    const metadata meta[] = {
        meta_builder().type(column_type::bigint).build(),
        meta_builder().type(column_type::char_).build(),
        meta_builder().type(column_type::blob).build(),
    };

    BOOST_TEST(map_metadata(map, 0, meta).type() == column_type::char_);
    BOOST_TEST(map_metadata(map, 1, meta).type() == column_type::bigint);
    BOOST_TEST(map_metadata(map, 2, meta).type() == column_type::blob);
}

BOOST_AUTO_TEST_CASE(map_field_view_)
{
    const std::array<std::size_t, 3> map{
        {1, 0, 2}
    };
    const auto fv = make_fv_arr(10, "abc", nullptr);

    BOOST_TEST(map_field_view(map, 0, fv) == field_view("abc"));
    BOOST_TEST(map_field_view(map, 1, fv) == field_view(10));
    BOOST_TEST(map_field_view(map, 2, fv) == field_view());
}

BOOST_AUTO_TEST_SUITE_END()
