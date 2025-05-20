//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_collection_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_integration/metadata_validator.hpp"

using namespace boost::mysql::test;

#define MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(getter) \
    {                                             \
        #getter, &boost::mysql::metadata::getter  \
    }

static struct flag_entry
{
    const char* name;
    meta_validator::flag_getter getter;
} flag_names[] = {
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_not_null),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_primary_key),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_unique_key),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_multiple_key),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_unsigned),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_zerofill),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_auto_increment),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(has_no_default_value),
    MYSQL_TEST_FLAG_GETTER_NAME_ENTRY(is_set_to_now_on_update)
};

BOOST_TEST_DONT_PRINT_LOG_VALUE(std::vector<flag_entry>::iterator)

static bool contains(
    meta_validator::flag_getter flag,
    const std::vector<meta_validator::flag_getter>& flagvec
)
{
    return std::find(flagvec.begin(), flagvec.end(), flag) != flagvec.end();
}

void meta_validator::validate(const boost::mysql::metadata& value) const
{
    // Fixed fields
    BOOST_TEST(value.database() == "boost_mysql_integtests");
    BOOST_TEST(value.table() == table_);
    BOOST_TEST(value.original_table() == org_table_);
    BOOST_TEST(value.column_name() == field_);
    BOOST_TEST(value.original_column_name() == org_field_);
    BOOST_TEST(value.column_length() > 0u);
    BOOST_TEST(value.type() == type_);
    BOOST_TEST(value.decimals() == decimals_);

    // Flags
    std::vector<flag_entry> all_flags(std::begin(flag_names), std::end(flag_names));

    for (flag_getter true_flag : flags_)
    {
        auto it = std::find_if(all_flags.begin(), all_flags.end(), [true_flag](const flag_entry& entry) {
            return entry.getter == true_flag;
        });
        BOOST_TEST_REQUIRE(it != all_flags.end());                // no repeated flag
        BOOST_TEST_REQUIRE(!contains(true_flag, ignore_flags_));  // ignore flags cannot be set to true
        BOOST_TEST((value.*true_flag)(), it->name);
        all_flags.erase(it);
    }

    for (const auto& entry : all_flags)
    {
        if (!contains(entry.getter, ignore_flags_))
        {
            BOOST_TEST(!(value.*entry.getter)(), entry.name);
        }
    }
}

void boost::mysql::test::validate_meta(
    const metadata_collection_view& actual,
    const std::vector<meta_validator>& expected
)
{
    BOOST_TEST_REQUIRE(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        expected[i].validate(actual[i]);
    }
}

void boost::mysql::test::validate_2fields_meta(metadata_collection_view fields, string_view table)
{
    validate_meta(
        fields,
        {meta_validator(table, "id", column_type::int_),
         meta_validator(table, "field_varchar", column_type::varchar)}
    );
}
