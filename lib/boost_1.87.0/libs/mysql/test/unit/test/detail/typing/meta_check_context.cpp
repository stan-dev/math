//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/typing/meta_check_context.hpp>
#include <boost/mysql/detail/typing/pos_map.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>

#include "test_common/printing.hpp"
#include "test_unit/create_meta.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::mysql::detail::meta_check_context;
using boost::mysql::detail::name_table_t;
using boost::mysql::detail::pos_absent;

BOOST_AUTO_TEST_SUITE(test_meta_check_context)

// Get the error message from the context
std::string get_errors(const meta_check_context& ctx)
{
    diagnostics diag;
    ctx.check_errors(diag);
    return diag.client_message();
}

BOOST_AUTO_TEST_CASE(column_type_to_str_)
{
    struct
    {
        column_type type;
        bool is_unsigned;
        string_view expected;
    } test_cases[] = {
        {column_type::tinyint,         false, "TINYINT"              },
        {column_type::tinyint,         true,  "TINYINT UNSIGNED"     },
        {column_type::smallint,        false, "SMALLINT"             },
        {column_type::smallint,        true,  "SMALLINT UNSIGNED"    },
        {column_type::mediumint,       false, "MEDIUMINT"            },
        {column_type::mediumint,       true,  "MEDIUMINT UNSIGNED"   },
        {column_type::int_,            false, "INT"                  },
        {column_type::int_,            true,  "INT UNSIGNED"         },
        {column_type::bigint,          false, "BIGINT"               },
        {column_type::bigint,          true,  "BIGINT UNSIGNED"      },
        {column_type::year,            false, "YEAR"                 },
        {column_type::year,            true,  "YEAR"                 },
        {column_type::float_,          false, "FLOAT"                },
        {column_type::float_,          true,  "FLOAT"                },
        {column_type::double_,         false, "DOUBLE"               },
        {column_type::double_,         true,  "DOUBLE"               },
        {column_type::date,            false, "DATE"                 },
        {column_type::datetime,        false, "DATETIME"             },
        {column_type::timestamp,       false, "TIMESTAMP"            },
        {column_type::time,            false, "TIME"                 },
        {column_type::char_,           false, "CHAR"                 },
        {column_type::varchar,         false, "VARCHAR"              },
        {column_type::text,            false, "TEXT"                 },
        {column_type::enum_,           false, "ENUM"                 },
        {column_type::set,             false, "SET"                  },
        {column_type::decimal,         false, "DECIMAL"              },
        {column_type::json,            false, "JSON"                 },
        {column_type::binary,          false, "BINARY"               },
        {column_type::varbinary,       false, "VARBINARY"            },
        {column_type::blob,            false, "BLOB"                 },
        {column_type::geometry,        false, "GEOMETRY"             },
        {static_cast<column_type>(76), true,  "<unknown column type>"}
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.expected)
        {
            auto meta = meta_builder().type(tc.type).unsigned_flag(tc.is_unsigned).build();
            BOOST_TEST(detail::column_type_to_str(meta) == tc.expected);
        }
    }
}

// Common data
const metadata meta[] = {
    meta_builder().type(column_type::bigint).build(),
    meta_builder().type(column_type::bit).build(),
    meta_builder().type(column_type::varchar).build(),
    meta_builder().type(column_type::blob).build(),
};
const string_view names[] = {"f1", "f2", "f3"};

BOOST_AUTO_TEST_CASE(accessors_fields_present)
{
    const std::size_t pos_map[] = {2, 0, 1};
    meta_check_context ctx(pos_map, names, meta);

    // Field 0
    BOOST_TEST(ctx.current_meta().type() == column_type::varchar);
    BOOST_TEST(!ctx.is_current_field_absent());
    ctx.advance();

    // Field 1
    BOOST_TEST(ctx.current_meta().type() == column_type::bigint);
    BOOST_TEST(!ctx.is_current_field_absent());
    ctx.advance();

    // Field 2
    BOOST_TEST(ctx.current_meta().type() == column_type::bit);
    BOOST_TEST(!ctx.is_current_field_absent());
}

BOOST_AUTO_TEST_CASE(accessors_some_fields_absent)
{
    const std::size_t pos_map[] = {1, pos_absent, 0};
    meta_check_context ctx(pos_map, name_table_t(), meta);

    // Field 0
    BOOST_TEST(ctx.current_meta().type() == column_type::bit);
    BOOST_TEST(!ctx.is_current_field_absent());
    ctx.advance();

    // Field 1
    BOOST_TEST(ctx.is_current_field_absent());
    ctx.advance();

    // Field 2
    BOOST_TEST(ctx.current_meta().type() == column_type::bigint);
    BOOST_TEST(!ctx.is_current_field_absent());
}

BOOST_AUTO_TEST_CASE(accessors_all_fields_absent)
{
    const std::size_t pos_map[] = {pos_absent, pos_absent, pos_absent};
    meta_check_context ctx(pos_map, name_table_t(), metadata_collection_view());

    // Field 0
    BOOST_TEST(ctx.is_current_field_absent());
    ctx.advance();

    // Field 1
    BOOST_TEST(ctx.is_current_field_absent());
    ctx.advance();

    // Field 2
    BOOST_TEST(ctx.is_current_field_absent());
}

BOOST_AUTO_TEST_CASE(nullability)
{
    const std::size_t pos_map[] = {0, 1, 2};
    meta_check_context ctx(pos_map, name_table_t(), meta);

    // Nullability not checked by default
    BOOST_TEST(!ctx.nullability_checked());

    // Explicitly setting it works
    ctx.set_nullability_checked();
    BOOST_TEST(ctx.nullability_checked());

    // Advancing resets it
    ctx.advance();
    BOOST_TEST(!ctx.nullability_checked());

    // Advancing again does nothing
    ctx.advance();
    BOOST_TEST(!ctx.nullability_checked());
}

BOOST_AUTO_TEST_CASE(add_field_absent_error_named)
{
    const std::size_t pos_map[] = {1, pos_absent, 0};
    const char* expected = "Field 'f2' is not present in the data returned by the server";

    meta_check_context ctx(pos_map, names, meta);
    ctx.advance();
    ctx.add_field_absent_error();
    BOOST_TEST(get_errors(ctx) == expected);
}

BOOST_AUTO_TEST_CASE(add_field_absent_error_unnamed)
{
    const std::size_t pos_map[] = {0, pos_absent};
    const char* expected =
        "Field in position 1 can't be mapped: there are more fields in your C++ data type than in your query";

    meta_check_context ctx(pos_map, name_table_t(), meta);
    ctx.advance();
    ctx.add_field_absent_error();
    BOOST_TEST(get_errors(ctx) == expected);
}

BOOST_AUTO_TEST_CASE(add_incompatible_types_error_named)
{
    const std::size_t pos_map[] = {1, 2, 0};
    const char* expected =
        "Incompatible types for field 'f2': C++ type 'cpp_type' is not compatible with DB type 'VARCHAR'";

    meta_check_context ctx(pos_map, names, meta);
    ctx.advance();
    ctx.add_type_mismatch_error("cpp_type");
    BOOST_TEST(get_errors(ctx) == expected);
}

BOOST_AUTO_TEST_CASE(add_incompatible_types_error_unnamed)
{
    const std::size_t pos_map[] = {0, 1, 2};
    const char* expected =
        "Incompatible types for field in position 1: C++ type 'other_type' is not compatible with DB type "
        "'BIT'";

    meta_check_context ctx(pos_map, name_table_t(), meta);
    ctx.advance();
    ctx.add_type_mismatch_error("other_type");
    BOOST_TEST(get_errors(ctx) == expected);
}

BOOST_AUTO_TEST_CASE(add_nullability_error_named)
{
    const std::size_t pos_map[] = {1, 2, 0};
    const char* expected =
        "NULL checks failed for field 'f1': the database type may be NULL, but the C++ type cannot. "
        "Use std::optional<T> or boost::optional<T>";

    meta_check_context ctx(pos_map, names, meta);
    ctx.add_nullability_error();
    BOOST_TEST(get_errors(ctx) == expected);
}

BOOST_AUTO_TEST_CASE(add_nullability_error_unnamed)
{
    const std::size_t pos_map[] = {0, 1, 2};
    const char* expected =
        "NULL checks failed for field in position 0: the database type may be NULL, but the C++ type "
        "cannot. Use std::optional<T> or boost::optional<T>";

    meta_check_context ctx(pos_map, name_table_t(), meta);
    ctx.add_nullability_error();
    BOOST_TEST(get_errors(ctx) == expected);
}

BOOST_AUTO_TEST_CASE(several_errors)
{
    const std::size_t pos_map[] = {3, pos_absent, 0};
    // clang-format off
    const char* expected = 
        "NULL checks failed for field 'f1': the database type may be NULL, but the C++ type cannot. Use std::optional<T> or boost::optional<T>\n"
        "Field 'f2' is not present in the data returned by the server\n"
        "Incompatible types for field 'f3': C++ type 'cpp_type' is not compatible with DB type 'BIGINT'\n"
        "NULL checks failed for field 'f3': the database type may be NULL, but the C++ type cannot. Use std::optional<T> or boost::optional<T>";
    // clang-format on

    meta_check_context ctx(pos_map, names, meta);
    ctx.add_nullability_error();
    ctx.advance();
    ctx.add_field_absent_error();
    ctx.advance();
    ctx.add_type_mismatch_error("cpp_type");
    ctx.add_nullability_error();
    BOOST_TEST(get_errors(ctx) == expected);
}

BOOST_AUTO_TEST_CASE(check_errors_no_error)
{
    const std::size_t pos_map[] = {0, 1};
    meta_check_context ctx(pos_map, name_table_t(), meta);
    diagnostics diag;
    error_code err = ctx.check_errors(diag);
    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.client_message() == "");
    BOOST_TEST(diag.server_message() == "");
}

BOOST_AUTO_TEST_CASE(check_errors_with_error)
{
    const std::size_t pos_map[] = {0, 1};
    const char* expected_msg =
        "Incompatible types for field in position 0: C++ type 'cpp_type' is not compatible with DB type "
        "'BIGINT'";

    meta_check_context ctx(pos_map, name_table_t(), meta);
    ctx.add_type_mismatch_error("cpp_type");
    diagnostics diag;
    error_code err = ctx.check_errors(diag);
    BOOST_TEST(err == client_errc::metadata_check_failed);
    BOOST_TEST(diag.client_message() == expected_msg);
    BOOST_TEST(diag.server_message() == "");
}

BOOST_AUTO_TEST_SUITE_END()
