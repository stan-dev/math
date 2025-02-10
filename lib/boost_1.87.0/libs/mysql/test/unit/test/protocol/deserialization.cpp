//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_categories.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/mysql_collations.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/protocol/capabilities.hpp>
#include <boost/mysql/impl/internal/protocol/deserialization.hpp>

#include <boost/test/unit_test.hpp>

#include <memory>

#include "operators.hpp"
#include "serialization_test.hpp"
#include "test_common/assert_buffer_equals.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/create_diagnostics.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_err.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_ok_frame.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/printing.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
namespace collations = boost::mysql::mysql_collations;
using boost::span;
using boost::mysql::client_errc;
using boost::mysql::column_type;
using boost::mysql::common_server_errc;
using boost::mysql::date;
using boost::mysql::datetime;
using boost::mysql::diagnostics;
using boost::mysql::error_code;
using boost::mysql::field_view;
using boost::mysql::get_mariadb_server_category;
using boost::mysql::get_mysql_server_category;
using boost::mysql::metadata;
using boost::mysql::string_view;

BOOST_TEST_DONT_PRINT_LOG_VALUE(execute_response::type_t)
BOOST_TEST_DONT_PRINT_LOG_VALUE(row_message::type_t)
BOOST_TEST_DONT_PRINT_LOG_VALUE(handhake_server_response::type_t)

BOOST_AUTO_TEST_SUITE(test_deserialization)

//
// OK packets
//
BOOST_AUTO_TEST_CASE(ok_view_success)
{
    struct
    {
        const char* name;
        ok_view expected;
        deserialization_buffer serialized;
    } test_cases[] = {
        // clang-format off
        {
            "successful_update",
            ok_builder()
                .affected_rows(4)
                .last_insert_id(0)
                .flags(34)
                .warnings(0)
                .info("Rows matched: 5  Changed: 4  Warnings: 0")
                .build(),
            {0x04, 0x00, 0x22, 0x00, 0x00, 0x00, 0x28, 0x52, 0x6f, 0x77, 0x73, 0x20, 0x6d, 0x61, 0x74, 0x63,
             0x68, 0x65, 0x64, 0x3a, 0x20, 0x35, 0x20, 0x20, 0x43, 0x68, 0x61, 0x6e, 0x67, 0x65, 0x64, 0x3a,
             0x20, 0x34, 0x20, 0x20, 0x57, 0x61, 0x72, 0x6e, 0x69, 0x6e, 0x67, 0x73, 0x3a, 0x20, 0x30},
        },
        {
            "successful_insert",
            ok_builder()
                .affected_rows(1)
                .last_insert_id(6)
                .flags(2)
                .warnings(0)
                .info("")
                .build(),
            {0x01, 0x06, 0x02, 0x00, 0x00, 0x00},
        },
        {
            "successful_login",
            ok_builder()
                .affected_rows(0)
                .last_insert_id(0)
                .flags(0x02)
                .warnings(0)
                .info("")
                .build(),
            {0x00, 0x00, 0x02, 0x00, 0x00, 0x00},
        }
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            ok_view actual{};
            error_code err = deserialize_ok_packet(tc.serialized, actual);

            // No error
            BOOST_TEST(err == error_code());

            // Actual value
            BOOST_TEST(actual.affected_rows == tc.expected.affected_rows);
            BOOST_TEST(actual.last_insert_id == tc.expected.last_insert_id);
            BOOST_TEST(actual.status_flags == tc.expected.status_flags);
            BOOST_TEST(actual.warnings == tc.expected.warnings);
            BOOST_TEST(actual.info == tc.expected.info);
        }
    }
}

BOOST_AUTO_TEST_CASE(ok_view_error)
{
    struct
    {
        const char* name;
        client_errc expected_err;
        deserialization_buffer serialized;
    } test_cases[] = {
        {"empty",                client_errc::incomplete_message, {}                                                    },
        {"error_affected_rows",  client_errc::incomplete_message, {0xff}                                                },
        {"error_last_insert_id", client_errc::incomplete_message, {0x01, 0xff}                                          },
        {"error_last_insert_id", client_errc::incomplete_message, {0x01, 0x06, 0x02}                                    },
        {"error_warnings",       client_errc::incomplete_message, {0x01, 0x06, 0x02, 0x00, 0x00}                        },
        {"error_info",           client_errc::incomplete_message, {0x04, 0x00, 0x22, 0x00, 0x00, 0x00, 0x28}            },
        {"extra_bytes",          client_errc::extra_bytes,        {0x01, 0x06, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00}}
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            ok_view value{};
            error_code err = deserialize_ok_packet(tc.serialized, value);
            BOOST_TEST(err == tc.expected_err);
        }
    }
}

//
// error packets
//
BOOST_AUTO_TEST_CASE(err_view_success)
{
    struct
    {
        const char* name;
        err_view expected;
        deserialization_buffer serialized;
        bool has_sql_state;
    } test_cases[] = {
        // clang-format off
        {
            "wrong_use_database",
            {1049, "Unknown database 'a'"},
            {0x19, 0x04, 0x23, 0x34, 0x32, 0x30, 0x30, 0x30, 0x55, 0x6e, 0x6b, 0x6e, 0x6f, 0x77,
             0x6e, 0x20, 0x64, 0x61, 0x74, 0x61, 0x62, 0x61, 0x73, 0x65, 0x20, 0x27, 0x61, 0x27},
            true,
        },
        {
            "unknown_table",
            {1146, "Table 'awesome.unknown' doesn't exist"},
            {0x7a, 0x04, 0x23, 0x34, 0x32, 0x53, 0x30, 0x32, 0x54, 0x61, 0x62, 0x6c, 0x65, 0x20, 0x27,
             0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x2e, 0x75, 0x6e, 0x6b, 0x6e, 0x6f, 0x77, 0x6e,
             0x27, 0x20, 0x64, 0x6f, 0x65, 0x73, 0x6e, 0x27, 0x74, 0x20, 0x65, 0x78, 0x69, 0x73, 0x74},
            true,
        },
        {
            "failed_login",
            {1045, "Access denied for user 'root'@'localhost' (using password: YES)"},
            {0x15, 0x04, 0x23, 0x32, 0x38, 0x30, 0x30, 0x30, 0x41, 0x63, 0x63, 0x65, 0x73, 0x73, 0x20,
             0x64, 0x65, 0x6e, 0x69, 0x65, 0x64, 0x20, 0x66, 0x6f, 0x72, 0x20, 0x75, 0x73, 0x65, 0x72,
             0x20, 0x27, 0x72, 0x6f, 0x6f, 0x74, 0x27, 0x40, 0x27, 0x6c, 0x6f, 0x63, 0x61, 0x6c, 0x68,
             0x6f, 0x73, 0x74, 0x27, 0x20, 0x28, 0x75, 0x73, 0x69, 0x6e, 0x67, 0x20, 0x70, 0x61, 0x73,
             0x73, 0x77, 0x6f, 0x72, 0x64, 0x3a, 0x20, 0x59, 0x45, 0x53, 0x29},
            true,
        },
        {
            "no_error_message",
            {1045, ""},
            {0x15, 0x04, 0x23, 0x32, 0x38, 0x30, 0x30, 0x30},
            true,
        },
        {
            "nosqlstate_too_many_connections",
            {1040, "Too many connections"},
            {0x10, 0x04, 0x54, 0x6f, 0x6f, 0x20, 0x6d, 0x61, 0x6e, 0x79, 0x20,
             0x63, 0x6f, 0x6e, 0x6e, 0x65, 0x63, 0x74, 0x69, 0x6f, 0x6e, 0x73},
            false,
        },
        {
            "nosqlstate_empty_err_msg",
            {1040, ""},
            {0x10, 0x04},
            false,
        },
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            err_view actual{};
            error_code err = deserialize_error_packet(tc.serialized, actual, tc.has_sql_state);

            // No error
            BOOST_TEST(err == error_code());

            // Actual value
            BOOST_TEST(actual.error_code == tc.expected.error_code);
            BOOST_TEST(actual.error_message == tc.expected.error_message);
        }
    }
}

BOOST_AUTO_TEST_CASE(err_view_error)
{
    struct
    {
        const char* name;
        deserialization_buffer serialized;
        bool has_sql_state;
    } test_cases[] = {
        {"empty",                       {},                       true },
        {"error_error_code",            {0x15},                   true },
        {"error_sql_state_marker",      {0x15, 0x04},             true },
        {"error_sql_state",             {0x15, 0x04, 0x23, 0x32}, true },
        {"nosqlstate_empty",            {},                       false},
        {"nosqlstate_error_error_code", {0x15},                   false},
    };
    // Note: not possible to get extra bytes here, since the last field is a string_eof

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            err_view value{};
            error_code err = deserialize_error_packet(tc.serialized, value, tc.has_sql_state);
            BOOST_TEST(err == client_errc::incomplete_message);
        }
    }
}

BOOST_AUTO_TEST_CASE(process_error_packet_)
{
    // It's OK to use err_builder() here, since the deserialization function
    // has already been tested
    struct
    {
        const char* name;
        db_flavor flavor;
        deserialization_buffer serialized;
        error_code ec;
        const char* msg;
    } test_cases[] = {
        {"bad_error_packet", db_flavor::mariadb, {0xff, 0x00, 0x01}, client_errc::incomplete_message, ""},
        {"code_lt_min",
         db_flavor::mariadb,
         err_builder().code(999).message("abc").build_body_without_header(),
         error_code(999, get_mariadb_server_category()),
         "abc"},
        {"code_common",
         db_flavor::mariadb,
         err_builder().code(1064).message("abc").build_body_without_header(),
         common_server_errc::er_parse_error,
         "abc"},
        {"code_common_hole_mysql",
         db_flavor::mysql,
         err_builder().code(1076).build_body_without_header(),
         error_code(1076, get_mysql_server_category()),
         ""},
        {"code_common_hole_mariadb",
         db_flavor::mariadb,
         err_builder().code(1076).build_body_without_header(),
         error_code(1076, get_mariadb_server_category()),
         ""},
        {"code_mysql",
         db_flavor::mysql,
         err_builder().code(4004).build_body_without_header(),
         error_code(4004, get_mysql_server_category()),
         ""},
        {"code_mariadb",
         db_flavor::mariadb,
         err_builder().code(4004).build_body_without_header(),
         error_code(4004, get_mariadb_server_category()),
         ""},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            diagnostics diag;
            auto ec = process_error_packet(tc.serialized, tc.flavor, diag);
            BOOST_TEST(ec == tc.ec);
            BOOST_TEST(diag.server_message() == tc.msg);
        }
    }
}

//
// coldef
//
BOOST_AUTO_TEST_CASE(coldef_view_success)
{
    struct
    {
        const char* name;
        coldef_view expected;
        deserialization_buffer serialized;
    } test_cases[] = {
        // clang-format off
        {
            "numeric_auto_increment_primary_key",
            meta_builder()
                .database("awesome")
                .table("test_table")
                .org_table("test_table")
                .name("id")
                .org_name("id")
                .collation_id(collations::binary)
                .column_length(11)
                .type(column_type::int_)
                .flags(
                    column_flags::not_null | column_flags::pri_key | column_flags::auto_increment |
                    column_flags::part_key
                )
                .decimals(0)
                .build_coldef(),
            {
                0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x0a, 0x74,
                0x65, 0x73, 0x74, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0a, 0x74, 0x65, 0x73, 0x74,
                0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x02, 0x69, 0x64, 0x02, 0x69, 0x64, 0x0c, 0x3f,
                0x00, 0x0b, 0x00, 0x00, 0x00, 0x03, 0x03, 0x42, 0x00, 0x00, 0x00
            },
        },
        {
            "varchar_field_aliased_field_and_table_names_join",
            meta_builder()
                .database("awesome")
                .table("child")
                .org_table("child_table")
                .name("field_alias")
                .org_name("field_varchar")
                .collation_id(collations::utf8_general_ci)
                .column_length(765)
                .type(column_type::varchar)
                .flags(0)
                .decimals(0)
                .build_coldef(),
            {
                0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
                0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68, 0x69,
                0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64,
                0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0b, 0x66,
                0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69,
                0x61, 0x73, 0x0d, 0x66, 0x69, 0x65, 0x6c, 0x64,
                0x5f, 0x76, 0x61, 0x72, 0x63, 0x68, 0x61, 0x72,
                0x0c, 0x21, 0x00, 0xfd, 0x02, 0x00, 0x00, 0xfd,
                0x00, 0x00, 0x00, 0x00, 0x00
            },
        },
        {
            "float_field",
            meta_builder()
                .database("awesome")
                .table("test_table")
                .org_table("test_table")
                .name("field_float")
                .org_name("field_float")
                .collation_id(collations::binary)
                .column_length(12)
                .type(column_type::float_)
                .flags(0)
                .decimals(31)
                .build_coldef(),
            {
                0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
                0x73, 0x6f, 0x6d, 0x65, 0x0a, 0x74, 0x65, 0x73,
                0x74, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0a,
                0x74, 0x65, 0x73, 0x74, 0x5f, 0x74, 0x61, 0x62,
                0x6c, 0x65, 0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64,
                0x5f, 0x66, 0x6c, 0x6f, 0x61, 0x74, 0x0b, 0x66,
                0x69, 0x65, 0x6c, 0x64, 0x5f, 0x66, 0x6c, 0x6f,
                0x61, 0x74, 0x0c, 0x3f, 0x00, 0x0c, 0x00, 0x00,
                0x00, 0x04, 0x00, 0x00, 0x1f, 0x00, 0x00
            },
        },
        {
            "no_final_padding", // edge case
            meta_builder()
                .database("awesome")
                .table("test_table")
                .org_table("test_table")
                .name("field_float")
                .org_name("field_float")
                .collation_id(collations::binary)
                .column_length(12)
                .type(column_type::float_)
                .flags(0)
                .decimals(31)
                .build_coldef(),
            {
                0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
                0x73, 0x6f, 0x6d, 0x65, 0x0a, 0x74, 0x65, 0x73,
                0x74, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0a,
                0x74, 0x65, 0x73, 0x74, 0x5f, 0x74, 0x61, 0x62,
                0x6c, 0x65, 0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64,
                0x5f, 0x66, 0x6c, 0x6f, 0x61, 0x74, 0x0b, 0x66,
                0x69, 0x65, 0x6c, 0x64, 0x5f, 0x66, 0x6c, 0x6f,
                0x61, 0x74, 0x0a, 0x3f, 0x00, 0x0c, 0x00, 0x00,
                0x00, 0x04, 0x00, 0x00, 0x1f
            },
        },
        {
            "more_final_padding", // test for extensibility - we don't fail if mysql adds more fields in the end
            meta_builder()
                .database("awesome")
                .table("test_table")
                .org_table("test_table")
                .name("field_float")
                .org_name("field_float")
                .collation_id(collations::binary)
                .column_length(12)
                .type(column_type::float_)
                .flags(0)
                .decimals(31)
                .build_coldef(),
            {
                0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
                0x73, 0x6f, 0x6d, 0x65, 0x0a, 0x74, 0x65, 0x73,
                0x74, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0a,
                0x74, 0x65, 0x73, 0x74, 0x5f, 0x74, 0x61, 0x62,
                0x6c, 0x65, 0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64,
                0x5f, 0x66, 0x6c, 0x6f, 0x61, 0x74, 0x0b, 0x66,
                0x69, 0x65, 0x6c, 0x64, 0x5f, 0x66, 0x6c, 0x6f,
                0x61, 0x74, 0x0d, 0x3f, 0x00, 0x0c, 0x00, 0x00,
                0x00, 0x04, 0x00, 0x00, 0x1f, 0x00, 0x00, 0x00
            },
        },
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            coldef_view actual{};
            error_code err = deserialize_column_definition(tc.serialized, actual);

            // No error
            BOOST_TEST_REQUIRE(err == error_code());

            // Actual value
            BOOST_TEST(actual.database == tc.expected.database);
            BOOST_TEST(actual.table == tc.expected.table);
            BOOST_TEST(actual.org_table == tc.expected.org_table);
            BOOST_TEST(actual.name == tc.expected.name);
            BOOST_TEST(actual.org_name == tc.expected.org_name);
            BOOST_TEST(actual.collation_id == tc.expected.collation_id);
            BOOST_TEST(actual.column_length == tc.expected.column_length);
            BOOST_TEST(actual.type == tc.expected.type);
            BOOST_TEST(actual.flags == tc.expected.flags);
            BOOST_TEST(actual.decimals == tc.expected.decimals);
        }
    }
}

BOOST_AUTO_TEST_CASE(coldef_view_error)
{
    struct
    {
        const char* name;
        error_code expected_err;
        deserialization_buffer serialized;
    } test_cases[] = {
        // clang-format off
        {
            "empty",
            client_errc::incomplete_message,
            {}
        },
        {
            "error_catalog",
            client_errc::incomplete_message,
            {0xff}
        },
        {
            "error_database",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0xff}
        },
        {
            "error_table",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0xff}
        },
        {   "error_org_table",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05,
            0x63, 0x68, 0x69, 0x6c, 0x64, 0xff}
        },
        {
            "error_name",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68, 0x69,
            0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0xff}
        },
        {
            "error_org_name",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68,
            0x69, 0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65,
            0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69, 0x61, 0x73, 0xff}
        },
        {
            "error_fixed_fields",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68,
            0x69, 0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65,
            0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69, 0x61, 0x73, 0x0d, 0x66, 0x69,
            0x65, 0x6c, 0x64, 0x5f, 0x76, 0x61, 0x72, 0x63, 0x68, 0x61, 0x72, 0xff}
        },
        {
            "error_collation_id",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68,
            0x69, 0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65,
            0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69, 0x61, 0x73, 0x0d, 0x66, 0x69,
            0x65, 0x6c, 0x64, 0x5f, 0x76, 0x61, 0x72, 0x63, 0x68, 0x61, 0x72, 0x01, 0x00}
        },
        {
            "error_column_length",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68,
            0x69, 0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65,
            0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69, 0x61, 0x73, 0x0d, 0x66, 0x69,
            0x65, 0x6c, 0x64, 0x5f, 0x76, 0x61, 0x72, 0x63, 0x68, 0x61, 0x72, 0x03, 0x00, 0x00, 0x00}
        },
        {
            "error_column_type",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63, 0x68, 0x69,
            0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0b, 0x66,
            0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69, 0x61, 0x73, 0x0d, 0x66, 0x69, 0x65, 0x6c, 0x64,
            0x5f, 0x76, 0x61, 0x72, 0x63, 0x68, 0x61, 0x72, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}
        },
        {
            "error_flags",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05,
            0x63, 0x68, 0x69, 0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74,
            0x61, 0x62, 0x6c, 0x65, 0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c,
            0x69, 0x61, 0x73, 0x0d, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x76, 0x61, 0x72,
            0x63, 0x68, 0x61, 0x72, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}
        },
        {
            "error_decimals",
            client_errc::incomplete_message,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x05, 0x63,
            0x68, 0x69, 0x6c, 0x64, 0x0b, 0x63, 0x68, 0x69, 0x6c, 0x64, 0x5f, 0x74, 0x61, 0x62,
            0x6c, 0x65, 0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x61, 0x6c, 0x69, 0x61, 0x73,
            0x0d, 0x66, 0x69, 0x65, 0x6c, 0x64, 0x5f, 0x76, 0x61, 0x72, 0x63, 0x68, 0x61, 0x72,
            0x09, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}
        },
        {
            "extra_bytes",
            client_errc::extra_bytes,
            {0x03, 0x64, 0x65, 0x66, 0x07, 0x61, 0x77, 0x65,
            0x73, 0x6f, 0x6d, 0x65, 0x0a, 0x74, 0x65, 0x73,
            0x74, 0x5f, 0x74, 0x61, 0x62, 0x6c, 0x65, 0x0a,
            0x74, 0x65, 0x73, 0x74, 0x5f, 0x74, 0x61, 0x62,
            0x6c, 0x65, 0x0b, 0x66, 0x69, 0x65, 0x6c, 0x64,
            0x5f, 0x66, 0x6c, 0x6f, 0x61, 0x74, 0x0b, 0x66,
            0x69, 0x65, 0x6c, 0x64, 0x5f, 0x66, 0x6c, 0x6f,
            0x61, 0x74, 0x0d, 0x3f, 0x00, 0x0c, 0x00, 0x00,
            0x00, 0x04, 0x00, 0x00, 0x1f, 0x00, 0x00, 0x00, 0xff}
        }
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            coldef_view value{};
            error_code err = deserialize_column_definition(tc.serialized, value);
            BOOST_TEST(err == tc.expected_err);
        }
    }
}

// OK response (ping & reset connection)
BOOST_AUTO_TEST_CASE(deserialize_ok_response_)
{
    struct
    {
        const char* name;
        deserialization_buffer message;
        error_code expected_err;
        const char* expected_msg;
        bool expected_backslash_escapes;
    } test_cases[] = {
        {"success",                      create_ok_body(ok_builder().build()),                error_code(),                      "", true },
        {"success_no_backslash_escapes",
         create_ok_body(ok_builder().no_backslash_escapes(true).build()),
         error_code(),
         "",                                                                                                                         false},
        {"empty_message",                {},                                                  client_errc::incomplete_message,   "", true },
        {"invalid_message_type",         {0xab},                                              client_errc::protocol_value_error, "", true },
        {"bad_ok_packet",                {0x00, 0x01},                                        client_errc::incomplete_message,   "", true },
        {"err_packet",
         err_builder().code(common_server_errc::er_bad_db_error).message("abc").build_body(),
         common_server_errc::er_bad_db_error,
         "abc",                                                                                                                      true },
        {"bad_err_packet",               {0xff, 0x01},                                        client_errc::incomplete_message,   "", true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            diagnostics diag;
            bool backslash_escapes = true;
            auto err = deserialize_ok_response(tc.message, db_flavor::mariadb, diag, backslash_escapes);

            BOOST_TEST(err == tc.expected_err);
            BOOST_TEST(diag.server_message() == tc.expected_msg);
            BOOST_TEST(backslash_escapes == tc.expected_backslash_escapes);
        }
    }
}

BOOST_AUTO_TEST_CASE(deserialize_prepare_stmt_response_impl_success)
{
    // Data (statement_id, num fields, num params)
    prepare_stmt_response expected{1, 2, 3};
    deserialization_buffer serialized{0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00};
    prepare_stmt_response actual{};
    auto err = deserialize_prepare_stmt_response_impl(serialized, actual);

    // No error
    BOOST_TEST_REQUIRE(err == error_code());

    // Actual value
    BOOST_TEST(actual.id == expected.id);
    BOOST_TEST(actual.num_columns == expected.num_columns);
    BOOST_TEST(actual.num_params == expected.num_params);
}

BOOST_AUTO_TEST_CASE(deserialize_prepare_stmt_response_impl_error)
{
    struct
    {
        const char* name;
        error_code expected_err;
        deserialization_buffer serialized;
    } test_cases[] = {
        {"empty",              client_errc::incomplete_message, {}                                              },
        {"error_id",           client_errc::incomplete_message, {0x01}                                          },
        {"error_num_columns",  client_errc::incomplete_message, {0x01, 0x00, 0x00, 0x00, 0x02}                  },
        {"error_num_params",   client_errc::incomplete_message, {0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03}      },
        {"error_reserved",     client_errc::incomplete_message, {0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x00}},
        {"error_num_warnings",
         client_errc::incomplete_message,
         {0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00}                                           },
        {"extra_bytes",
         client_errc::extra_bytes,
         {0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0xff}                               },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            prepare_stmt_response output{};
            auto err = deserialize_prepare_stmt_response_impl(tc.serialized, output);
            BOOST_TEST(err == tc.expected_err);
        }
    }
}

BOOST_AUTO_TEST_CASE(deserialize_prepare_stmt_response_success)
{
    // Data (statement_id, num fields, num params)
    prepare_stmt_response expected{1, 2, 3};
    deserialization_buffer serialized{0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00};
    prepare_stmt_response actual{};
    diagnostics diag;

    auto err = deserialize_prepare_stmt_response(serialized, db_flavor::mysql, actual, diag);

    // No error
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(diag == diagnostics());

    // Actual value
    BOOST_TEST(actual.id == expected.id);
    BOOST_TEST(actual.num_columns == expected.num_columns);
    BOOST_TEST(actual.num_params == expected.num_params);
}

BOOST_AUTO_TEST_CASE(deserialize_prepare_stmt_response_error)
{
    struct
    {
        const char* name;
        error_code expected_err;
        const char* expected_diag;
        deserialization_buffer serialized;
    } test_cases[] = {
        // clang-format off
        {
            "error_message_type",
            client_errc::incomplete_message,
            "",
            {},
        },
        {
            "unknown_message_type",
            client_errc::protocol_value_error,
            "",
            {0xab, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00},
        },
        {
            "error_packet",
            common_server_errc::er_bad_db_error,
            "bad db",
            err_builder().code(common_server_errc::er_bad_db_error).message("bad db").build_body(),
        },
        {
            "error_deserializing_response",
            client_errc::incomplete_message,
            "",
            {0x00, 0x01, 0x00},
        },
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            prepare_stmt_response output{};
            diagnostics diag;

            auto err = deserialize_prepare_stmt_response(tc.serialized, db_flavor::mariadb, output, diag);

            BOOST_TEST(err == tc.expected_err);
            BOOST_TEST(diag.server_message() == tc.expected_diag);
        }
    }
}

//
// execute response
//
BOOST_AUTO_TEST_CASE(deserialize_execute_response_ok_packet)
{
    deserialization_buffer serialized{0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00};
    diagnostics diag;

    auto response = deserialize_execute_response(serialized, db_flavor::mariadb, diag);

    BOOST_TEST_REQUIRE(response.type == execute_response::type_t::ok_packet);
    BOOST_TEST(response.data.ok_pack.affected_rows == 0u);
    BOOST_TEST(response.data.ok_pack.status_flags == 2u);
}

BOOST_AUTO_TEST_CASE(deserialize_execute_response_num_fields)
{
    struct
    {
        const char* name;
        deserialization_buffer serialized;
        std::size_t num_fields;
    } test_cases[] = {
        {"1",                    {0x01},             1     },
        {"0xfa",                 {0xfa},             0xfa  },
        {"0xfb_no_local_infile", {0xfb},             0xfb  }, // legal when LOCAL INFILE capability not enabled
        {"0xfb_local_infile",    {0xfc, 0xfb, 0x00}, 0xfb  }, // sent LOCAL INFILE capability is enabled
        {"0xff",                 {0xfc, 0xff, 0x00}, 0xff  },
        {"0x01ff",               {0xfc, 0xff, 0x01}, 0x01ff},
        {"max",                  {0xfc, 0xff, 0xff}, 0xffff},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            diagnostics diag;

            auto response = deserialize_execute_response(tc.serialized, db_flavor::mysql, diag);

            BOOST_TEST_REQUIRE(response.type == execute_response::type_t::num_fields);
            BOOST_TEST(response.data.num_fields == tc.num_fields);
            BOOST_TEST(diag.server_message() == "");
        }
    }
}

BOOST_AUTO_TEST_CASE(deserialize_execute_response_error)
{
    struct
    {
        const char* name;
        deserialization_buffer serialized;
        error_code err;
        const char* expected_info;
    } test_cases[] = {
        {"server_error",
         {0xff, 0x7a, 0x04, 0x23, 0x34, 0x32, 0x53, 0x30, 0x32, 0x54, 0x61, 0x62, 0x6c, 0x65,
          0x20, 0x27, 0x6d, 0x79, 0x74, 0x65, 0x73, 0x74, 0x2e, 0x61, 0x62, 0x63, 0x27, 0x20,
          0x64, 0x6f, 0x65, 0x73, 0x6e, 0x27, 0x74, 0x20, 0x65, 0x78, 0x69, 0x73, 0x74},
         common_server_errc::er_no_such_table,
         "Table 'mytest.abc' doesn't exist"                                                                                   },
        {"bad_server_error", {0xff, 0x00},                                               client_errc::incomplete_message,   ""},
        {"bad_ok_packet",    {0x00, 0xff},                                               client_errc::incomplete_message,   ""},
        {"bad_num_fields",   {0xfc, 0xff, 0x00, 0x01},                                   client_errc::extra_bytes,          ""},
        {"zero_num_fields",  {0xfc, 0x00, 0x00},                                         client_errc::protocol_value_error, ""},
        {"3byte_integer",    {0xfd, 0xff, 0xff, 0xff},                                   client_errc::protocol_value_error, ""},
        {"8byte_integer",
         {0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
         client_errc::protocol_value_error,
         ""                                                                                                                   },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            diagnostics diag;

            auto response = deserialize_execute_response(tc.serialized, db_flavor::mysql, diag);

            BOOST_TEST_REQUIRE(response.type == execute_response::type_t::error);
            BOOST_TEST(response.data.err == tc.err);
            BOOST_TEST(diag.server_message() == tc.expected_info);
        }
    }
}

//
// row message
//
BOOST_AUTO_TEST_CASE(deserialize_row_message_row)
{
    deserialization_buffer serialized{create_text_row_body("abc", 10)};
    diagnostics diag;

    auto response = deserialize_row_message(serialized, db_flavor::mysql, diag);

    BOOST_TEST_REQUIRE(response.type == row_message::type_t::row);
    BOOST_TEST(response.data.row.data() == serialized.data());
    BOOST_TEST(response.data.row.size() == serialized.size());
}

BOOST_AUTO_TEST_CASE(deserialize_row_message_ok_packet)
{
    deserialization_buffer serialized{
        create_eof_body(ok_builder().affected_rows(42).last_insert_id(1).info("abc").build())
    };
    diagnostics diag;

    auto response = deserialize_row_message(serialized, db_flavor::mysql, diag);

    BOOST_TEST_REQUIRE(response.type == row_message::type_t::ok_packet);
    BOOST_TEST(response.data.ok_pack.affected_rows == 42u);
    BOOST_TEST(response.data.ok_pack.last_insert_id == 1u);
    BOOST_TEST(response.data.ok_pack.info == "abc");
}

BOOST_AUTO_TEST_CASE(deserialize_row_message_error)
{
    struct
    {
        const char* name;
        deserialization_buffer serialized;
        error_code expected_error;
        const char* expected_info;
    } test_cases[] = {
        // clang-format off
        {
            "invalid_ok_packet",
            { 0xfe, 0x00, 0x00, 0x02, 0x00, 0x00 }, // 1 byte missing
            client_errc::incomplete_message,
            ""
        },
        {
            "error_packet",
            {
                0xff, 0x19, 0x04, 0x23, 0x34, 0x32, 0x30, 0x30, 0x30, 0x55, 0x6e, 0x6b,
                0x6e, 0x6f, 0x77, 0x6e, 0x20, 0x64, 0x61, 0x74,
                0x61, 0x62, 0x61, 0x73, 0x65, 0x20, 0x27, 0x61, 0x27
            },
            common_server_errc::er_bad_db_error,
            "Unknown database 'a'"
        },
        {
            "invalid_error_packet",
            { 0xff, 0x19 }, // bytes missing
            client_errc::incomplete_message,
            ""
        },
        {
            "empty_message",
            {},
            client_errc::incomplete_message,
            ""
        }
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            diagnostics diag;

            auto msg = deserialize_row_message(tc.serialized, db_flavor::mysql, diag);

            BOOST_TEST_REQUIRE(msg.type == row_message::type_t::error);
            BOOST_TEST(msg.data.err == tc.expected_error);
            BOOST_TEST(diag.server_message() == tc.expected_info);
        }
    }
}

//
// deserialize_row
//
BOOST_AUTO_TEST_CASE(deserialize_row_success)
{
    struct
    {
        const char* name;
        resultset_encoding encoding;
        deserialization_buffer serialized;
        std::vector<field_view> expected;
        std::vector<metadata> meta;
    } test_cases[] = {
        // clang-format off
        // Text
        {
            "text_one_value",
            resultset_encoding::text,
            {0x01, 0x35},
            make_fv_vector(std::int64_t(5)),
            create_metas({ column_type::tinyint })
        },
        {
            "text_one_null",
            resultset_encoding::text,
            {0xfb},
            make_fv_vector(nullptr),
            create_metas({ column_type::tinyint} )
        },
        {
            "text_several_values",
            resultset_encoding::text,
            {0x03, 0x76, 0x61, 0x6c, 0x02, 0x32, 0x31, 0x03, 0x30, 0x2e, 0x30},
            make_fv_vector("val", std::int64_t(21), 0.0f),
            create_metas({ column_type::varchar, column_type::int_, column_type::float_})
        },
        {
            "text_several_values_one_null",
            resultset_encoding::text,
            {0x03, 0x76, 0x61, 0x6c, 0xfb, 0x03, 0x76, 0x61, 0x6c},
            make_fv_vector("val", nullptr, "val"),
            create_metas({ column_type::varchar, column_type::int_, column_type::varchar })
        },
        {
            "text_several_nulls",
            resultset_encoding::text,
            {0xfb, 0xfb, 0xfb},
            make_fv_vector(nullptr, nullptr, nullptr),
            create_metas({ column_type::varchar, column_type::int_, column_type::datetime })
        },

        // Binary
        {
            "binary_one_value",
            resultset_encoding::binary,
            {0x00, 0x00, 0x14},
            make_fv_vector(std::int64_t(20)),
            create_metas({ column_type::tinyint })
        },
        {
            "binary_one_null",
            resultset_encoding::binary,
            {0x00, 0x04},
            make_fv_vector(nullptr),
            create_metas({ column_type::tinyint })
        },
        {
            "binary_two_values",
            resultset_encoding::binary,
            {0x00, 0x00, 0x03, 0x6d, 0x69, 0x6e, 0x6d, 0x07},
            make_fv_vector("min", std::int64_t(1901)),
            create_metas({ column_type::varchar, column_type::smallint })
        },
        {
            "binary_one_value_one_null",
            resultset_encoding::binary,
            {0x00, 0x08, 0x03, 0x6d, 0x61, 0x78},
            make_fv_vector("max", nullptr),
            create_metas({ column_type::varchar, column_type::tinyint })
        },
        {
            "binary_two_nulls",
            resultset_encoding::binary,
            {0x00, 0x0c},
            make_fv_vector(nullptr, nullptr),
            create_metas({ column_type::tinyint, column_type::tinyint })
        },
        {
            "binary_six_nulls",
            resultset_encoding::binary,
            {0x00, 0xfc},
            std::vector<field_view>(6, field_view()),
            std::vector<metadata>(6, create_meta(column_type::tinyint))
        },
        {
            "binary_seven_nulls",
            resultset_encoding::binary,
            {0x00, 0xfc, 0x01},
            std::vector<field_view>(7, field_view()),
            std::vector<metadata>(7, create_meta(column_type::tinyint))
        },
        {
            "binary_several_values",
            resultset_encoding::binary,
            {
                0x00, 0x90, 0x00, 0xfd, 0x03, 0x61, 0x62, 0x63,
                0xc3, 0xf5, 0x48, 0x40, 0x02, 0x61, 0x62, 0x04,
                0xe2, 0x07, 0x0a, 0x05, 0x71, 0x99, 0x6d, 0xe2,
                0x93, 0x4d, 0xf5, 0x3d
            },
            make_fv_vector(
                std::int64_t(-3),
                "abc",
                nullptr,
                3.14f,
                "ab",
                nullptr,
                date(2018u, 10u, 5u),
                3.10e-10
            ),
            create_metas({
                column_type::tinyint,
                column_type::varchar,
                column_type::int_,
                column_type::float_,
                column_type::char_,
                column_type::int_,
                column_type::date,
                column_type::double_,
            })
        }
    };
    // clang-format on

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Allocate exactly what is expected, to facilitate tooling for overrun detection
            std::unique_ptr<field_view[]> actual{new field_view[tc.expected.size()]};
            span<field_view> actual_span{actual.get(), tc.expected.size()};

            auto err = deserialize_row(tc.encoding, tc.serialized, tc.meta, actual_span);

            BOOST_TEST_REQUIRE(err == error_code());
            std::vector<field_view> actual_vec{actual_span.begin(), actual_span.end()};
            BOOST_TEST(actual_vec == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(deserialize_row_error)
{
    struct
    {
        const char* name;
        resultset_encoding encoding;
        deserialization_buffer serialized;
        error_code expected;
        std::vector<metadata> meta;
    } test_cases[] = {
        // clang-format off
        // text
        {
            "text_no_space_string_single",
            resultset_encoding::text,
            {0x02, 0x00},
            client_errc::incomplete_message,
            create_metas({ column_type::smallint })
        },
        {
            "text_no_space_string_final",
            resultset_encoding::text,
            {0x01, 0x35, 0x02, 0x35},
            client_errc::incomplete_message,
            create_metas({ column_type::tinyint, column_type::smallint }),
        },
        {
            "text_no_space_null_single",
            resultset_encoding::text,
            {},
            client_errc::incomplete_message,
            create_metas({ column_type::tinyint })
        },
        {
            "text_no_space_null_final",
            resultset_encoding::text,
            {0xfb},
            client_errc::incomplete_message,
            create_metas({ column_type::tinyint, column_type::tinyint }),
        },
        {
            "text_extra_bytes",
            resultset_encoding::text,
            {0x01, 0x35, 0xfb, 0x00},
            client_errc::extra_bytes,
            create_metas({ column_type::tinyint, column_type::tinyint })
        },
        {
            "text_contained_value_error_single",
            resultset_encoding::text,
            {0x01, 0x00},
            client_errc::protocol_value_error,
            create_metas({ column_type::date })
        },
        {
            "text_contained_value_error_middle",
            resultset_encoding::text,
            {0xfb, 0x01, 0x00, 0xfb},
            client_errc::protocol_value_error,
            create_metas({ column_type::date, column_type::date, column_type::date })
        },
        {
            "text_row_for_empty_meta",
            resultset_encoding::text,
            {0xfb, 0x01, 0x00, 0xfb},
            client_errc::extra_bytes,
            {}
        },

        // binary
        {
            "binary_no_space_null_bitmap_1",
            resultset_encoding::binary,
            {0x00},
            client_errc::incomplete_message,
            create_metas({ column_type::tinyint })
        },
        {
            "binary_no_space_null_bitmap_2",
            resultset_encoding::binary,
            {0x00, 0xfc},
            client_errc::incomplete_message,
            std::vector<metadata>(7, create_meta(column_type::tinyint))
        },
        {
            "binary_no_space_value_single",
            resultset_encoding::binary,
            {0x00, 0x00},
            client_errc::incomplete_message,
            create_metas({ column_type::tinyint })
        },
        {
            "binary_no_space_value_last",
            resultset_encoding::binary,
            {0x00, 0x00, 0x01},
            client_errc::incomplete_message,
            create_metas({ column_type::tinyint, column_type::tinyint })
        },
        {
            "binary_no_space_value_middle",
            resultset_encoding::binary,
            {0x00, 0x00, 0x01},
            client_errc::incomplete_message,
            create_metas({ column_type::tinyint, column_type::tinyint, column_type::tinyint })
        },
        {
            "binary_extra_bytes",
            resultset_encoding::binary,
            {0x00, 0x00, 0x01, 0x02},
            client_errc::extra_bytes,
            create_metas({ column_type::tinyint })
        },
        {
            "binary_row_for_empty_meta",
            resultset_encoding::binary,
            {0xfb, 0x01, 0x00, 0xfb},
            client_errc::extra_bytes,
            {}
        },
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Allocate exactly what is expected, to facilitate tooling for overrun detection
            std::unique_ptr<field_view[]> actual{new field_view[tc.meta.size()]};
            span<field_view> actual_span{actual.get(), tc.meta.size()};

            auto err = deserialize_row(tc.encoding, tc.serialized, tc.meta, actual_span);

            BOOST_TEST(err == tc.expected);
        }
    }
}

//
// server hello
//
BOOST_AUTO_TEST_CASE(deserialize_server_hello_impl_success)
{
    // Data
    constexpr std::uint8_t auth_plugin_data[] = {0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f, 0x03, 0x5a,
                                                 0x74, 0x05, 0x28, 0x2b, 0x7f, 0x21, 0x43, 0x4a, 0x21, 0x62};

    constexpr std::uint32_t caps = CLIENT_LONG_PASSWORD | CLIENT_FOUND_ROWS | CLIENT_LONG_FLAG |
                                   CLIENT_CONNECT_WITH_DB | CLIENT_NO_SCHEMA | CLIENT_COMPRESS | CLIENT_ODBC |
                                   CLIENT_LOCAL_FILES | CLIENT_IGNORE_SPACE | CLIENT_PROTOCOL_41 |
                                   CLIENT_INTERACTIVE | CLIENT_IGNORE_SIGPIPE | CLIENT_TRANSACTIONS |
                                   CLIENT_RESERVED |           // old flag, but set in this frame
                                   CLIENT_SECURE_CONNECTION |  // old flag, but set in this frame
                                   CLIENT_MULTI_STATEMENTS | CLIENT_MULTI_RESULTS | CLIENT_PS_MULTI_RESULTS |
                                   CLIENT_PLUGIN_AUTH | CLIENT_CONNECT_ATTRS |
                                   CLIENT_PLUGIN_AUTH_LENENC_CLIENT_DATA |
                                   CLIENT_CAN_HANDLE_EXPIRED_PASSWORDS | CLIENT_SESSION_TRACK |
                                   CLIENT_DEPRECATE_EOF | CLIENT_REMEMBER_OPTIONS;

    deserialization_buffer serialized{0x35, 0x2e, 0x37, 0x2e, 0x32, 0x37, 0x2d, 0x30, 0x75, 0x62, 0x75, 0x6e,
                                      0x74, 0x75, 0x30, 0x2e, 0x31, 0x39, 0x2e, 0x30, 0x34, 0x2e, 0x31, 0x00,
                                      0x02, 0x00, 0x00, 0x00, 0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
                                      0x00, 0xff, 0xf7, 0x08, 0x02, 0x00, 0xff, 0x81, 0x15, 0x00, 0x00, 0x00,
                                      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x5a, 0x74, 0x05, 0x28,
                                      0x2b, 0x7f, 0x21, 0x43, 0x4a, 0x21, 0x62, 0x00, 0x6d, 0x79, 0x73, 0x71,
                                      0x6c, 0x5f, 0x6e, 0x61, 0x74, 0x69, 0x76, 0x65, 0x5f, 0x70, 0x61, 0x73,
                                      0x73, 0x77, 0x6f, 0x72, 0x64, 0x00};

    // Deserialize
    server_hello actual{};
    auto err = deserialize_server_hello_impl(serialized, actual);

    // No error
    BOOST_TEST_REQUIRE(err == error_code());

    // Actual value
    BOOST_TEST(actual.server == db_flavor::mysql);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(actual.auth_plugin_data.to_span(), auth_plugin_data);
    BOOST_TEST(actual.server_capabilities == capabilities(caps));
    BOOST_TEST(actual.auth_plugin_name == "mysql_native_password");

    // TODO: mysql8, mariadb, edge case where auth plugin length is < 13
}

BOOST_AUTO_TEST_CASE(deserialize_server_hello_impl_error)
{
    struct
    {
        const char* name;
        deserialization_buffer serialized;
        error_code expected_err;
    } test_cases[] = {
        // clang-format off
        {
            "error_server_version",
            {0x10, 0x11}, // no NULL-terminator
            client_errc::incomplete_message
        },  
        {
            "error_connection_id",
            {0x2e, 0x31, 0x00,
            0x02},
            client_errc::incomplete_message
        },
        {
            "error_auth_plugin_data_1",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a},
            client_errc::incomplete_message
        },
        {
            "error_auth_plugin_data_filler",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f},
            client_errc::incomplete_message
        },
        {
            "error_capability_flags_low",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff},
            client_errc::incomplete_message
        },
        {
            "error_character_set",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7},
            client_errc::incomplete_message
        },
        {
            "error_status_flags",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7,
            0x08,
            0x02},
            client_errc::incomplete_message
        },
        {
            "error_capability_flags_high",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7,
            0x08,
            0x02, 0x00,
            0xff},
            client_errc::incomplete_message
        },
        {
            "error_auth_plugin_data_length",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7,
            0x08,
            0x02, 0x00,
            0xff, 0x81},
            client_errc::incomplete_message
        },
        {
            "error_reserved",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7,
            0x08,
            0x02, 0x00,
            0xff, 0x81,
            0x15,
            0x00, 0x00},
            client_errc::incomplete_message
        },
        {
            "error_auth_plugin_data_2",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7,
            0x08, // charset
            0x02, 0x00,
            0xff, 0x81, // caps high
            0x15, // auth plugin data len
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // reserved
            // plugin data 2 has 1 byte less
            0x03, 0x5a, 0x74, 0x05, 0x28, 0x2b, 0x7f, 0x21, 0x43, 0x4a, 0x21, 0x62},
            client_errc::incomplete_message
        },
        {
            "error_auth_plugin_name",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7,
            0x08, // charset
            0x02, 0x00,
            0xff, 0x81, // caps high
            0x15, // auth plugin data len
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // reserved
            0x03, 0x5a, 0x74, 0x05, 0x28, 0x2b, 0x7f, 0x21, 0x43, 0x4a, 0x21, 0x62, 0x00, // plugin data 2
            // plugin name
            0x6d, 0x79, 0x73, 0x71, 0x6c, 0x5f, 0x6e, 0x61, 0x74,
            0x69, 0x76, 0x65, 0x5f, 0x70, 0x61, 0x73, 0x73, 0x77, 0x6f, 0x72, 0x64}, // no NULL terminator
            client_errc::incomplete_message
        },
        {
            "no_plugin_auth_capability",
            {0x2e, 0x31, 0x00,
            0x02, 0x00, 0x00, 0x00,
            0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f,
            0x00,
            0xff, 0xf7,
            0x08, // charset
            0x02, 0x00,
            0xf7, 0x81, // caps high - no CLIENT_PLUGIN_AUTH
            0x00, // auth plugin data len
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // reserved
            0x03, 0x5a, 0x74, 0x05, 0x28, 0x2b, 0x7f, 0x21, 0x43, 0x4a, 0x21, 0x62, 0x00, // plugin data 2
            },
            client_errc::server_unsupported
        }
        // clang-format on
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            server_hello value{};
            auto err = deserialize_server_hello_impl(tc.serialized, value);
            BOOST_TEST(err == tc.expected_err);
        }
    }
}

//
// deserialize_server_hello
//
BOOST_AUTO_TEST_CASE(deserialize_server_hello_success)
{
    // Data
    constexpr std::uint8_t auth_plugin_data[] = {0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70, 0x2f, 0x03, 0x5a,
                                                 0x74, 0x05, 0x28, 0x2b, 0x7f, 0x21, 0x43, 0x4a, 0x21, 0x62};

    constexpr std::uint32_t caps = CLIENT_LONG_PASSWORD | CLIENT_FOUND_ROWS | CLIENT_LONG_FLAG |
                                   CLIENT_CONNECT_WITH_DB | CLIENT_NO_SCHEMA | CLIENT_COMPRESS | CLIENT_ODBC |
                                   CLIENT_LOCAL_FILES | CLIENT_IGNORE_SPACE | CLIENT_PROTOCOL_41 |
                                   CLIENT_INTERACTIVE | CLIENT_IGNORE_SIGPIPE | CLIENT_TRANSACTIONS |
                                   CLIENT_RESERVED |           // old flag, but set in this frame
                                   CLIENT_SECURE_CONNECTION |  // old flag, but set in this frame
                                   CLIENT_MULTI_STATEMENTS | CLIENT_MULTI_RESULTS | CLIENT_PS_MULTI_RESULTS |
                                   CLIENT_PLUGIN_AUTH | CLIENT_CONNECT_ATTRS |
                                   CLIENT_PLUGIN_AUTH_LENENC_CLIENT_DATA |
                                   CLIENT_CAN_HANDLE_EXPIRED_PASSWORDS | CLIENT_SESSION_TRACK |
                                   CLIENT_DEPRECATE_EOF | CLIENT_REMEMBER_OPTIONS;

    deserialization_buffer serialized{0x0a, 0x35, 0x2e, 0x37, 0x2e, 0x32, 0x37, 0x2d, 0x30, 0x75, 0x62, 0x75,
                                      0x6e, 0x74, 0x75, 0x30, 0x2e, 0x31, 0x39, 0x2e, 0x30, 0x34, 0x2e, 0x31,
                                      0x00, 0x02, 0x00, 0x00, 0x00, 0x52, 0x1a, 0x50, 0x3a, 0x4b, 0x12, 0x70,
                                      0x2f, 0x00, 0xff, 0xf7, 0x08, 0x02, 0x00, 0xff, 0x81, 0x15, 0x00, 0x00,
                                      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x5a, 0x74, 0x05,
                                      0x28, 0x2b, 0x7f, 0x21, 0x43, 0x4a, 0x21, 0x62, 0x00, 0x6d, 0x79, 0x73,
                                      0x71, 0x6c, 0x5f, 0x6e, 0x61, 0x74, 0x69, 0x76, 0x65, 0x5f, 0x70, 0x61,
                                      0x73, 0x73, 0x77, 0x6f, 0x72, 0x64, 0x00};

    // Deserialize
    server_hello actual{};
    diagnostics diag;
    auto err = deserialize_server_hello(serialized, actual, diag);

    // No error
    BOOST_TEST_REQUIRE(err == error_code());
    BOOST_TEST(diag == diagnostics());

    // Actual value
    BOOST_TEST(actual.server == db_flavor::mysql);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(actual.auth_plugin_data.to_span(), auth_plugin_data);
    BOOST_TEST(actual.server_capabilities == capabilities(caps));
    BOOST_TEST(actual.auth_plugin_name == "mysql_native_password");
}

BOOST_AUTO_TEST_CASE(deserialize_server_hello_error)
{
    struct
    {
        const char* name;
        deserialization_buffer serialized;
        error_code expected_err;
        diagnostics expected_diag;
    } test_cases[] = {
        // Regression check: we were failing to deserialize this correctly because of
        // missing sql_state
        {"err_packet",
         {0xff, 0x10, 0x04, 0x54, 0x6f, 0x6f, 0x20, 0x6d, 0x61, 0x6e, 0x79, 0x20,
          0x63, 0x6f, 0x6e, 0x6e, 0x65, 0x63, 0x74, 0x69, 0x6f, 0x6e, 0x73},
         common_server_errc::er_con_count_error,
         create_server_diag("Too many connections")},
        {"handshake_v9",
         {0x09, 0x6d, 0x79, 0x73, 0x71, 0x6c, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0a, 0x00},
         client_errc::server_unsupported,
         {}},
        {"bad_server_hello", {0x0a, 0x11}, client_errc::incomplete_message, {}},
        {"unknown_header", {0xab, 0x00, 0x01}, client_errc::protocol_value_error, {}},
        {"empty_message", {}, client_errc::incomplete_message, {}},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            server_hello msg{};
            diagnostics diag;

            // Deserialize
            auto err = deserialize_server_hello(tc.serialized, msg, diag);

            // Check
            BOOST_TEST(err == tc.expected_err);
            BOOST_TEST(diag == tc.expected_diag);
        }
    }
}

//
// auth switch
//
BOOST_AUTO_TEST_CASE(deserialize_auth_switch_success)
{
    // Data
    constexpr std::uint8_t auth_data[] = {0x49, 0x49, 0x7e, 0x51, 0x5d, 0x1f, 0x19, 0x6a, 0x0f, 0x5a,
                                          0x63, 0x15, 0x3e, 0x28, 0x31, 0x3e, 0x3c, 0x79, 0x09, 0x7c};

    deserialization_buffer serialized{0x6d, 0x79, 0x73, 0x71, 0x6c, 0x5f, 0x6e, 0x61, 0x74, 0x69, 0x76,
                                      0x65, 0x5f, 0x70, 0x61, 0x73, 0x73, 0x77, 0x6f, 0x72, 0x64, 0x00,
                                      0x49, 0x49, 0x7e, 0x51, 0x5d, 0x1f, 0x19, 0x6a, 0x0f, 0x5a, 0x63,
                                      0x15, 0x3e, 0x28, 0x31, 0x3e, 0x3c, 0x79, 0x09, 0x7c, 0x00};

    // Deserialize
    auth_switch actual{};
    auto err = deserialize_auth_switch(serialized, actual);

    // No error
    BOOST_TEST_REQUIRE(err == error_code());

    // Actual value
    BOOST_TEST(actual.plugin_name == "mysql_native_password");
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(actual.auth_data, auth_data);

    // TODO: error in plugin_name, extra bytes (?)
    // TODO: edge case where plugin data is empty or doesn't end in a 0
}

//
// handshake server response
//
BOOST_AUTO_TEST_CASE(deserialize_handshake_server_response_more_data)
{
    constexpr std::uint8_t auth_data[] = {0x61, 0x62, 0x63};

    deserialization_buffer serialized{0x01, 0x61, 0x62, 0x63};

    // Deserialize
    diagnostics diag;
    auto response = deserialize_handshake_server_response(serialized, db_flavor::mysql, diag);

    // Actual value
    BOOST_TEST_REQUIRE(response.type == handhake_server_response::type_t::auth_more_data);
    BOOST_MYSQL_ASSERT_BUFFER_EQUALS(response.data.more_data, auth_data);
}
// TODO: ok packet
// TODO: ok follows
// TODO: error packet
// TODO: auth switch
// TODO: error in message type, unknown message type, bad OK packet, bad auth switch

BOOST_AUTO_TEST_SUITE_END()
