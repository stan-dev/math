//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/coldef_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_unit/create_meta.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
namespace collations = boost::mysql::mysql_collations;
namespace column_flags = boost::mysql::detail::column_flags;

// TODO: these tests should be more targetted towards metadata's object responsibilities
// and less integration-y

BOOST_AUTO_TEST_SUITE(test_metadata)

BOOST_AUTO_TEST_CASE(int_primary_key)
{
    detail::coldef_view msg{
        "awesome",
        "test_table",
        "test_table",
        "id",
        "id",
        collations::binary,
        11,
        column_type::int_,
        column_flags::pri_key | column_flags::auto_increment | column_flags::not_null,
        0,
    };
    auto meta = detail::access::construct<metadata>(msg, true);

    BOOST_TEST(meta.database() == "awesome");
    BOOST_TEST(meta.table() == "test_table");
    BOOST_TEST(meta.original_table() == "test_table");
    BOOST_TEST(meta.column_name() == "id");
    BOOST_TEST(meta.original_column_name() == "id");
    BOOST_TEST(meta.column_length() == 11u);
    BOOST_TEST(meta.type() == column_type::int_);
    BOOST_TEST(meta.decimals() == 0u);
    BOOST_TEST(meta.is_not_null());
    BOOST_TEST(meta.is_primary_key());
    BOOST_TEST(!meta.is_unique_key());
    BOOST_TEST(!meta.is_multiple_key());
    BOOST_TEST(!meta.is_unsigned());
    BOOST_TEST(!meta.is_zerofill());
    BOOST_TEST(meta.is_auto_increment());
    BOOST_TEST(!meta.has_no_default_value());
    BOOST_TEST(!meta.is_set_to_now_on_update());
}

BOOST_AUTO_TEST_CASE(varchar_with_alias)
{
    detail::coldef_view msg{
        "awesome",
        "child",
        "child_table",
        "field_alias",
        "field_varchar",
        collations::utf8mb4_general_ci,
        765,
        column_type::varchar,
        0,
        0,
    };
    auto meta = detail::access::construct<metadata>(msg, true);

    BOOST_TEST(meta.database() == "awesome");
    BOOST_TEST(meta.table() == "child");
    BOOST_TEST(meta.original_table() == "child_table");
    BOOST_TEST(meta.column_name() == "field_alias");
    BOOST_TEST(meta.original_column_name() == "field_varchar");
    BOOST_TEST(meta.column_length() == 765u);
    BOOST_TEST(meta.type() == column_type::varchar);
    BOOST_TEST(meta.decimals() == 0u);
    BOOST_TEST(!meta.is_not_null());
    BOOST_TEST(!meta.is_primary_key());
    BOOST_TEST(!meta.is_unique_key());
    BOOST_TEST(!meta.is_multiple_key());
    BOOST_TEST(!meta.is_unsigned());
    BOOST_TEST(!meta.is_zerofill());
    BOOST_TEST(!meta.is_auto_increment());
    BOOST_TEST(!meta.has_no_default_value());
    BOOST_TEST(!meta.is_set_to_now_on_update());
}

BOOST_AUTO_TEST_CASE(float_)
{
    detail::coldef_view msg{
        "awesome",
        "test_table",
        "test_table",
        "field_float",
        "field_float",
        collations::binary,
        12,
        column_type::float_,
        0,
        31,
    };
    auto meta = detail::access::construct<metadata>(msg, true);

    BOOST_TEST(meta.database() == "awesome");
    BOOST_TEST(meta.table() == "test_table");
    BOOST_TEST(meta.original_table() == "test_table");
    BOOST_TEST(meta.column_name() == "field_float");
    BOOST_TEST(meta.original_column_name() == "field_float");
    BOOST_TEST(meta.column_length() == 12u);
    BOOST_TEST(meta.type() == column_type::float_);
    BOOST_TEST(meta.decimals() == 31u);
    BOOST_TEST(!meta.is_not_null());
    BOOST_TEST(!meta.is_primary_key());
    BOOST_TEST(!meta.is_unique_key());
    BOOST_TEST(!meta.is_multiple_key());
    BOOST_TEST(!meta.is_unsigned());
    BOOST_TEST(!meta.is_zerofill());
    BOOST_TEST(!meta.is_auto_increment());
    BOOST_TEST(!meta.has_no_default_value());
    BOOST_TEST(!meta.is_set_to_now_on_update());
}

BOOST_AUTO_TEST_CASE(dont_copy_strings)
{
    detail::coldef_view msg{
        "awesome",
        "child",
        "child_table",
        "field_alias",
        "field_varchar",
        collations::utf8mb4_general_ci,
        765,
        column_type::varchar,
        0,
        0,
    };
    auto meta = detail::access::construct<metadata>(msg, false);

    BOOST_TEST(meta.database() == "");
    BOOST_TEST(meta.table() == "");
    BOOST_TEST(meta.original_table() == "");
    BOOST_TEST(meta.column_name() == "");
    BOOST_TEST(meta.original_column_name() == "");
    BOOST_TEST(meta.column_length() == 765u);
    BOOST_TEST(meta.type() == column_type::varchar);
    BOOST_TEST(meta.decimals() == 0u);
    BOOST_TEST(!meta.is_not_null());
    BOOST_TEST(!meta.is_primary_key());
    BOOST_TEST(!meta.is_unique_key());
    BOOST_TEST(!meta.is_multiple_key());
    BOOST_TEST(!meta.is_unsigned());
    BOOST_TEST(!meta.is_zerofill());
    BOOST_TEST(!meta.is_auto_increment());
    BOOST_TEST(!meta.has_no_default_value());
    BOOST_TEST(!meta.is_set_to_now_on_update());
}

BOOST_AUTO_TEST_CASE(string_ownership)
{
    // Create the meta object
    std::string colname = "col1";
    auto msg = meta_builder().name(colname).build_coldef();
    auto meta = detail::access::construct<metadata>(msg, true);

    // Check that we actually copy the data
    colname = "abcd";
    BOOST_TEST(meta.column_name() == "col1");

    // TODO: the other strings
}

BOOST_AUTO_TEST_SUITE_END()  // test_metadata
