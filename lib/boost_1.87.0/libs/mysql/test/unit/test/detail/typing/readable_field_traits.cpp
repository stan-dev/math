//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/field.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/time.hpp>

#include <boost/mysql/detail/typing/meta_check_context.hpp>
#include <boost/mysql/detail/typing/pos_map.hpp>
#include <boost/mysql/detail/typing/readable_field_traits.hpp>
#include <boost/mysql/detail/typing/row_traits.hpp>

#include <boost/core/span.hpp>
#include <boost/optional/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <limits>
#include <string>

#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/custom_allocator.hpp"

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
#include <optional>
#endif

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::span;
using boost::mysql::detail::is_readable_field;
using boost::mysql::detail::meta_check_context;
using boost::mysql::detail::name_table_t;
using boost::mysql::detail::pos_absent;
using boost::mysql::detail::readable_field_traits;

// Don't attempt to print std::chrono values
BOOST_TEST_DONT_PRINT_LOG_VALUE(boost::mysql::time)

BOOST_AUTO_TEST_SUITE(test_readable_field_traits)

//
// readable_field
//

struct other_traits : std::char_traits<char>
{
};
using string_with_alloc = std::basic_string<char, std::char_traits<char>, custom_allocator<char>>;
using string_with_traits = std::basic_string<char, other_traits>;
using blob_with_alloc = std::vector<unsigned char, custom_allocator<unsigned char>>;
using detail::meta_check_field;

struct unrelated
{
};

static_assert(is_readable_field<unsigned char>::value, "");
static_assert(is_readable_field<signed char>::value, "");
static_assert(is_readable_field<char>::value, "");
static_assert(is_readable_field<short>::value, "");
static_assert(is_readable_field<unsigned short>::value, "");
static_assert(is_readable_field<int>::value, "");
static_assert(is_readable_field<unsigned int>::value, "");
static_assert(is_readable_field<long>::value, "");
static_assert(is_readable_field<unsigned long>::value, "");
static_assert(is_readable_field<long long>::value, "");
static_assert(is_readable_field<unsigned long long>::value, "");
static_assert(is_readable_field<std::int8_t>::value, "");
static_assert(is_readable_field<std::uint8_t>::value, "");
static_assert(is_readable_field<std::int16_t>::value, "");
static_assert(is_readable_field<std::uint16_t>::value, "");
static_assert(is_readable_field<std::int32_t>::value, "");
static_assert(is_readable_field<std::uint32_t>::value, "");
static_assert(is_readable_field<std::int64_t>::value, "");
static_assert(is_readable_field<std::uint64_t>::value, "");
static_assert(is_readable_field<std::int_least8_t>::value, "");
static_assert(is_readable_field<std::uint_least8_t>::value, "");
static_assert(is_readable_field<std::int_least16_t>::value, "");
static_assert(is_readable_field<std::uint_least16_t>::value, "");
static_assert(is_readable_field<std::int_least32_t>::value, "");
static_assert(is_readable_field<std::uint_least32_t>::value, "");
static_assert(is_readable_field<std::int_least64_t>::value, "");
static_assert(is_readable_field<std::uint_least64_t>::value, "");
static_assert(is_readable_field<std::int_fast8_t>::value, "");
static_assert(is_readable_field<std::uint_fast8_t>::value, "");
static_assert(is_readable_field<std::int_fast16_t>::value, "");
static_assert(is_readable_field<std::uint_fast16_t>::value, "");
static_assert(is_readable_field<std::int_fast32_t>::value, "");
static_assert(is_readable_field<std::uint_fast32_t>::value, "");
static_assert(is_readable_field<std::int_fast64_t>::value, "");
static_assert(is_readable_field<std::uint_fast64_t>::value, "");
static_assert(is_readable_field<float>::value, "");
static_assert(is_readable_field<double>::value, "");
static_assert(is_readable_field<boost::mysql::date>::value, "");
static_assert(is_readable_field<boost::mysql::datetime>::value, "");
static_assert(is_readable_field<boost::mysql::time>::value, "");
static_assert(is_readable_field<bool>::value, "");

// string types
static_assert(is_readable_field<std::string>::value, "");
static_assert(is_readable_field<string_with_alloc>::value, "");
static_assert(!is_readable_field<string_with_traits>::value, "");
static_assert(!is_readable_field<std::wstring>::value, "");
static_assert(!is_readable_field<string_view>::value, "");

// blob types
static_assert(is_readable_field<blob>::value, "");
static_assert(is_readable_field<blob_with_alloc>::value, "");
static_assert(!is_readable_field<blob_view>::value, "");

// references not accepted
static_assert(!is_readable_field<int&>::value, "");
static_assert(!is_readable_field<const int&>::value, "");
static_assert(!is_readable_field<int&&>::value, "");
static_assert(!is_readable_field<const int&>::value, "");
static_assert(!is_readable_field<std::string&>::value, "");

// optionals
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
static_assert(is_readable_field<std::optional<int>>::value, "");
static_assert(is_readable_field<std::optional<std::string>>::value, "");
#endif
static_assert(is_readable_field<boost::optional<blob>>::value, "");
static_assert(is_readable_field<boost::optional<datetime>>::value, "");
static_assert(!is_readable_field<boost::optional<void*>>::value, "");
static_assert(!is_readable_field<boost::optional<unrelated>>::value, "");

// other types not accepted
static_assert(!is_readable_field<std::nullptr_t>::value, "");
static_assert(!is_readable_field<field_view>::value, "");
static_assert(!is_readable_field<field>::value, "");
static_assert(!is_readable_field<const char*>::value, "");
static_assert(!is_readable_field<void*>::value, "");
static_assert(!is_readable_field<unrelated>::value, "");
static_assert(!is_readable_field<const field_view*>::value, "");

// const-qualified objects are not accepted
static_assert(!is_readable_field<const int>::value, "");
static_assert(!is_readable_field<const std::string>::value, "");
static_assert(!is_readable_field<boost::optional<const int>>::value, "");

BOOST_AUTO_TEST_SUITE(meta_check_field_)

using single_field_check_fn = void (*)(meta_check_context&);

constexpr struct cpp_type_descriptor
{
    const char* name;
    single_field_check_fn check_fn;
} cpp_type_descriptors[] = {
    {"int8_t",   &meta_check_field<std::int8_t>       },
    {"uint8_t",  &meta_check_field<std::uint8_t>      },
    {"int16_t",  &meta_check_field<std::int16_t>      },
    {"uint16_t", &meta_check_field<std::uint16_t>     },
    {"int32_t",  &meta_check_field<std::int32_t>      },
    {"uint32_t", &meta_check_field<std::uint32_t>     },
    {"int64_t",  &meta_check_field<std::int64_t>      },
    {"uint64_t", &meta_check_field<std::uint64_t>     },
    {"bool",     &meta_check_field<bool>              },
    {"float",    &meta_check_field<float>             },
    {"double",   &meta_check_field<double>            },
    {"date",     &meta_check_field<date>              },
    {"datetime", &meta_check_field<datetime>          },
    {"time",     &meta_check_field<boost::mysql::time>},
    {"string",   &meta_check_field<std::string>       },
    {"blob",     &meta_check_field<blob>              },
};
constexpr auto cpp_type_descriptors_size = sizeof(cpp_type_descriptors) / sizeof(cpp_type_descriptor);

constexpr struct db_type_descriptor
{
    const char* name;
    const char* pretty_name;
    column_type type;
    bool is_unsigned;
} db_type_descriptors[] = {
    {"TINYINT",               "tinyint",    column_type::tinyint,   false},
    {"TINYINT UNSIGNED",      "tinyintu",   column_type::tinyint,   true },
    {"SMALLINT",              "smallint",   column_type::smallint,  false},
    {"SMALLINT UNSIGNED",     "smallintu",  column_type::smallint,  true },
    {"MEDIUMINT",             "mediumint",  column_type::mediumint, false},
    {"MEDIUMINT UNSIGNED",    "mediumintu", column_type::mediumint, true },
    {"INT",                   "int",        column_type::int_,      false},
    {"INT UNSIGNED",          "intu",       column_type::int_,      true },
    {"BIGINT",                "bigint",     column_type::bigint,    false},
    {"BIGINT UNSIGNED",       "bigintu",    column_type::bigint,    true },
    {"YEAR",                  "year",       column_type::year,      true },
    {"BIT",                   "bit",        column_type::bit,       true },
    {"FLOAT",                 "float",      column_type::float_,    false},
    {"DOUBLE",                "double",     column_type::double_,   false},
    {"DATE",                  "date",       column_type::date,      false},
    {"DATETIME",              "datetime",   column_type::datetime,  false},
    {"TIMESTAMP",             "timestamp",  column_type::timestamp, false},
    {"TIME",                  "time",       column_type::time,      false},
    {"CHAR",                  "char",       column_type::char_,     false},
    {"VARCHAR",               "varchar",    column_type::varchar,   false},
    {"TEXT",                  "text",       column_type::text,      false},
    {"ENUM",                  "enum",       column_type::enum_,     false},
    {"SET",                   "set",        column_type::set,       false},
    {"JSON",                  "json",       column_type::json,      false},
    {"DECIMAL",               "decimal",    column_type::decimal,   false},
    {"BINARY",                "binary",     column_type::binary,    false},
    {"VARBINARY",             "varbinary",  column_type::varbinary, false},
    {"BLOB",                  "blob",       column_type::blob,      false},
    {"GEOMETRY",              "geometry",   column_type::geometry,  false},
    {"<unknown column type>", "unknown",    column_type::unknown,   false},
};
constexpr auto db_type_descriptors_size = sizeof(db_type_descriptors) / sizeof(db_type_descriptor);

struct compat_matrix_row
{
    std::array<bool, cpp_type_descriptors_size> compat;
};

// Looks like clang-format crashes when it sees this
// clang-format off
constexpr std::array<compat_matrix_row, db_type_descriptors_size> compat_matrix{{
    {{{1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0}}},  // TINYINT
    {{{0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}}},  // TINYINT UNSIGNED
    {{{0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},  // SMALLINT
    {{{0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}}},  // SMALLINT UNSIGNED
    {{{0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},  // MEDIUMINT
    {{{0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}}},  // MEDIUMINT UNSIGNED
    {{{0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},  // INT
    {{{0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}}},  // INT UNSIGNED
    {{{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},  // BIGINT
    {{{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}}},  // BIGINT UNSIGNED
    {{{0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}}},  // YEAR
    {{{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}}},  // BIT
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0}}},  // FLOAT
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}}},  // DOUBLE
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}}},  // DATE
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}}},  // DATETIME
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}}},  // TIMESTAMP
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}}},  // TIME
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}}},  // CHAR
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}}},  // VARCHAR
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}}},  // TEXT
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}}},  // ENUM
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}}},  // SET
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}}},  // JSON
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}}},  // DECIMAL
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}},  // BINARY
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}},  // VARBINARY
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}},  // BLOB
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}},  // GEOMETRY
    {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}},  // UNKNOWN
}};
// clang-format on

BOOST_AUTO_TEST_CASE(basic_types_compatible)
{
    for (std::size_t i = 0; i < db_type_descriptors_size; ++i)
    {
        const auto& row = compat_matrix.at(i);
        for (std::size_t j = 0; j < cpp_type_descriptors_size; ++j)
        {
            if (!row.compat.at(j))
                continue;
            auto db_desc = db_type_descriptors[i];
            auto cpp_desc = cpp_type_descriptors[j];
            BOOST_TEST_CONTEXT(db_desc.pretty_name << "_" << cpp_desc.name)
            {
                const std::size_t pos_map[] = {0};
                const metadata meta[] = {
                    meta_builder()
                        .type(db_desc.type)
                        .unsigned_flag(db_desc.is_unsigned)
                        .nullable(false)
                        .build(),
                };
                meta_check_context ctx(pos_map, name_table_t(), meta);

                cpp_desc.check_fn(ctx);
                diagnostics diag;
                auto err = ctx.check_errors(diag);
                BOOST_TEST(err == error_code());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(basic_types_incompatible)
{
    for (std::size_t i = 0; i < db_type_descriptors_size; ++i)
    {
        const auto& row = compat_matrix.at(i);
        for (std::size_t j = 0; j < cpp_type_descriptors_size; ++j)
        {
            if (row.compat.at(j))
                continue;
            auto db_desc = db_type_descriptors[i];
            auto cpp_desc = cpp_type_descriptors[j];

            BOOST_TEST_CONTEXT(db_desc.pretty_name << "_" << cpp_desc.name)
            {
                const std::size_t pos_map[] = {0};
                const metadata meta[] = {
                    meta_builder()
                        .type(db_desc.type)
                        .unsigned_flag(db_desc.is_unsigned)
                        .nullable(false)
                        .build(),
                };
                meta_check_context ctx(pos_map, name_table_t(), meta);

                cpp_desc.check_fn(ctx);
                diagnostics diag;
                auto err = ctx.check_errors(diag);

                BOOST_TEST(err == client_errc::metadata_check_failed);
                std::string msg = "Incompatible types for field in position 0: C++ type '" +
                                  std::string(cpp_desc.name) + "' is not compatible with DB type '" +
                                  std::string(db_desc.name) + "'";
                BOOST_TEST(diag.client_message() == msg);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(nullable_error)
{
    const metadata meta[] = {
        meta_builder().type(column_type::float_).unsigned_flag(false).nullable(true).build(),
    };
    const std::size_t pos_map[] = {0};
    meta_check_context ctx(pos_map, name_table_t(), meta);

    meta_check_field<double>(ctx);
    diagnostics diag;
    auto err = ctx.check_errors(diag);
    BOOST_TEST(err == client_errc::metadata_check_failed);
    BOOST_TEST(
        diag.client_message() ==
        "NULL checks failed for field in position 0: the database type may be NULL, but the C++ type cannot. "
        "Use std::optional<T> or boost::optional<T>"
    );
}

BOOST_AUTO_TEST_CASE(optionals)
{
    struct
    {
        const char* name;
        single_field_check_fn check_fn;
        bool nullable;
    } test_cases[] = {
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
        {"std_optional_not_nullable",   &meta_check_field<std::optional<double>>,   false},
        {"std_optional_nullable",       &meta_check_field<std::optional<double>>,   true },
#endif
        {"boost_optional_not_nullable", &meta_check_field<boost::optional<double>>, false},
        {"boost_optional_nullable",     &meta_check_field<boost::optional<double>>, true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            const std::size_t pos_map[] = {0};
            const metadata meta[] = {
                meta_builder().type(column_type::float_).unsigned_flag(false).nullable(tc.nullable).build(),
            };
            meta_check_context ctx(pos_map, name_table_t(), meta);

            tc.check_fn(ctx);
            diagnostics diag;
            auto err = ctx.check_errors(diag);
            BOOST_TEST(err == error_code());
        }
    }
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(parse_)

template <class T>
struct parse_and_check
{
    T expected;

    parse_and_check(T expected) : expected(std::move(expected)) {}

    void operator()(field_view from) const
    {
        T actual{};
        auto err = readable_field_traits<T>::parse(from, actual);
        BOOST_TEST(err == error_code());
        BOOST_TEST(actual == expected);
    }
};

BOOST_AUTO_TEST_CASE(success)
{
    struct
    {
        const char* name;
        field_view from;
        std::function<void(field_view)> fn;
    } test_cases[] = {
        {"int8_signed_regular", field_view(42), parse_and_check<std::int8_t>(42)},
        {"int8_signed_min", field_view(-0x80), parse_and_check<std::int8_t>(-0x80)},
        {"int8_signed_max", field_view(0x7f), parse_and_check<std::int8_t>(0x7f)},
        {"int8_unsigned_regular", field_view(42u), parse_and_check<std::int8_t>(42u)},
        {"int8_unsigned_max", field_view(0x7fu), parse_and_check<std::int8_t>(0x7fu)},

        {"uint8_regular", field_view(42u), parse_and_check<std::uint8_t>(42u)},
        {"uint8_min", field_view(0u), parse_and_check<std::uint8_t>(0u)},
        {"uint8_max", field_view(0xffu), parse_and_check<std::uint8_t>(0xffu)},

        {"int16_signed_regular", field_view(42), parse_and_check<std::int16_t>(42)},
        {"int16_signed_min", field_view(-0x8000), parse_and_check<std::int16_t>(-0x8000)},
        {"int16_signed_max", field_view(0x7f00), parse_and_check<std::int16_t>(0x7f00)},
        {"int16_unsigned_regular", field_view(42u), parse_and_check<std::int16_t>(42u)},
        {"int16_unsigned_max", field_view(0x7f00u), parse_and_check<std::int16_t>(0x7f00u)},

        {"uint16_regular", field_view(42u), parse_and_check<std::uint16_t>(42u)},
        {"uint16_min", field_view(0u), parse_and_check<std::uint16_t>(0u)},
        {"uint16_max", field_view(0xffffu), parse_and_check<std::uint16_t>(0xffffu)},

        {"int32_signed_regular", field_view(42), parse_and_check<std::int32_t>(42)},
        {"int32_signed_min", field_view(-0x80000000LL), parse_and_check<std::int32_t>(-0x80000000LL)},
        {"int32_signed_max", field_view(0x7f000000), parse_and_check<std::int32_t>(0x7f000000)},
        {"int32_unsigned_regular", field_view(42u), parse_and_check<std::int32_t>(42u)},
        {"int32_unsigned_max", field_view(0x7f000000u), parse_and_check<std::int32_t>(0x7f000000u)},

        {"uint32_regular", field_view(42u), parse_and_check<std::uint32_t>(42u)},
        {"uint32_min", field_view(0u), parse_and_check<std::uint32_t>(0u)},
        {"uint32_max", field_view(0xffffffffu), parse_and_check<std::uint32_t>(0xffffffffu)},

        {"int64_signed_regular", field_view(42), parse_and_check<std::int64_t>(42)},
        {"int64_signed_min",
         field_view(-0x7fffffffffffffff - 1),
         parse_and_check<std::int64_t>(-0x7fffffffffffffff - 1)},
        {"int64_signed_max", field_view(0x7f00000000000000), parse_and_check<std::int64_t>(0x7f00000000000000)
        },
        {"int64_unsigned_regular", field_view(42u), parse_and_check<std::int64_t>(42u)},
        {"int64_unsigned_max",
         field_view(0x7f00000000000000u),
         parse_and_check<std::int64_t>(0x7f00000000000000u)},

        {"uint64_regular", field_view(42u), parse_and_check<std::uint64_t>(42u)},
        {"uint64_min", field_view(0u), parse_and_check<std::uint64_t>(0u)},
        {"uint64_max", field_view(0xffffffffffffffffu), parse_and_check<std::uint64_t>(0xffffffffffffffffu)},

        {"bool_zero", field_view(0), parse_and_check<bool>(false)},
        {"bool_one", field_view(1), parse_and_check<bool>(true)},
        {"bool_other", field_view(2), parse_and_check<bool>(true)},

        {"float", field_view(4.2f), parse_and_check<float>(4.2f)},

        {"double_float", field_view(4.2f), parse_and_check<double>(4.2f)},
        {"double_double", field_view(4.2), parse_and_check<double>(4.2)},

        {"date", field_view(date(2020, 1, 2)), parse_and_check<date>(date(2020, 1, 2))},
        {"datetime", field_view(datetime(2020, 1, 2)), parse_and_check<datetime>(datetime(2020, 1, 2))},
        {"time", field_view(maket(10, 1, 1)), parse_and_check<boost::mysql::time>(maket(10, 1, 1))},
        {"string", field_view("abc"), parse_and_check<std::string>("abc")},
        {"blob", field_view(makebv("\0\1")), parse_and_check<blob>({0, 1})},

        {"boost_optional_empty",
         field_view(),
         parse_and_check<boost::optional<float>>(boost::optional<float>())},
        {"boost_optional_nonempty",
         field_view(4.2f),
         parse_and_check<boost::optional<float>>(boost::optional<float>(4.2f))},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name) { tc.fn(tc.from); }
    }
}

// std::optional, which doesn't implement equality or stream operators
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
BOOST_AUTO_TEST_CASE(std_optional_null)
{
    std::optional<int> actual = 42;
    auto err = readable_field_traits<std::optional<int>>::parse(field_view(), actual);
    BOOST_TEST(err == error_code());
    BOOST_TEST(!actual.has_value());
}

BOOST_AUTO_TEST_CASE(std_optional_not_null)
{
    std::optional<int> actual;
    auto err = readable_field_traits<std::optional<int>>::parse(field_view(42), actual);
    BOOST_TEST(err == error_code());
    BOOST_TEST_REQUIRE(actual.has_value());
    BOOST_TEST(actual.value() == 42);
}
#endif

template <class T>
error_code parse_and_discard(field_view f)
{
    T t;
    return readable_field_traits<T>::parse(f, t);
}

BOOST_AUTO_TEST_CASE(errors)
{
    field_view i64_absmin((std::numeric_limits<std::int64_t>::min)());
    field_view i64_absmax((std::numeric_limits<std::int64_t>::max)());
    field_view ui64_absmax((std::numeric_limits<std::uint64_t>::max)());

    struct
    {
        const char* name;
        field_view from;
        error_code (*parse_fn)(field_view);
    } test_cases[] = {
        {"int8_null",                       field_view(),                    parse_and_discard<std::int8_t>                 },
        {"int8_badtype",                    field_view(4.2),                 parse_and_discard<std::int8_t>                 },
        {"int8_signed_ltmin",               field_view(-0x81),               parse_and_discard<std::int8_t>                 },
        {"int8_signed_gtmax",               field_view(0x80),                parse_and_discard<std::int8_t>                 },
        {"int8_signed_absmin",              i64_absmin,                      parse_and_discard<std::int8_t>                 },
        {"int8_signed_absmax",              i64_absmax,                      parse_and_discard<std::int8_t>                 },
        {"int8_unsigned_gtmax",             field_view(0x80u),               parse_and_discard<std::int8_t>                 },
        {"int8_unsigned_absmax",            ui64_absmax,                     parse_and_discard<std::int8_t>                 },

        {"uint8_null",                      field_view(),                    parse_and_discard<std::uint8_t>                },
        {"uint8_badtype",                   field_view(4.2),                 parse_and_discard<std::uint8_t>                },
        {"uint8_gtmax",                     field_view(0x100u),              parse_and_discard<std::uint8_t>                },
        {"uint8_absmax",                    ui64_absmax,                     parse_and_discard<std::uint8_t>                },

        {"int16_null",                      field_view(),                    parse_and_discard<std::int16_t>                },
        {"int16_badtype",                   field_view(4.2),                 parse_and_discard<std::int16_t>                },
        {"int16_signed_ltmin",              field_view(-0x8001),             parse_and_discard<std::int16_t>                },
        {"int16_signed_gtmax",              field_view(0x8000),              parse_and_discard<std::int16_t>                },
        {"int16_signed_absmin",             i64_absmin,                      parse_and_discard<std::int16_t>                },
        {"int16_signed_absmax",             i64_absmax,                      parse_and_discard<std::int16_t>                },
        {"int16_unsigned_gtmax",            field_view(0x8000u),             parse_and_discard<std::int16_t>                },
        {"int16_unsigned_absmax",           ui64_absmax,                     parse_and_discard<std::int16_t>                },

        {"uint16_null",                     field_view(),                    parse_and_discard<std::uint16_t>               },
        {"uint16_badtype",                  field_view(4.2),                 parse_and_discard<std::uint16_t>               },
        {"uint16_gtmax",                    field_view(0x10000u),            parse_and_discard<std::uint16_t>               },
        {"uint16_absmax",                   ui64_absmax,                     parse_and_discard<std::uint16_t>               },

        {"int32_null",                      field_view(),                    parse_and_discard<std::int32_t>                },
        {"int32_badtype",                   field_view(4.2),                 parse_and_discard<std::int32_t>                },
        {"int32_signed_ltmin",              field_view(-0x80000001LL),       parse_and_discard<std::int32_t>                },
        {"int32_signed_gtmax",              field_view(0x80000000L),         parse_and_discard<std::int32_t>                },
        {"int32_signed_absmin",             i64_absmin,                      parse_and_discard<std::int32_t>                },
        {"int32_signed_absmax",             i64_absmax,                      parse_and_discard<std::int32_t>                },
        {"int32_unsigned_gtmax",            field_view(0x80000000uL),        parse_and_discard<std::int32_t>                },
        {"int32_unsigned_absmax",           ui64_absmax,                     parse_and_discard<std::int32_t>                },

        {"uint32_null",                     field_view(),                    parse_and_discard<std::uint32_t>               },
        {"uint32_badtype",                  field_view("abc"),               parse_and_discard<std::uint32_t>               },
        {"uint32_gtmax",                    field_view(0x100000000u),        parse_and_discard<std::uint32_t>               },
        {"uint32_absmax",                   ui64_absmax,                     parse_and_discard<std::uint32_t>               },

        {"int64_null",                      field_view(),                    parse_and_discard<std::int64_t>                },
        {"int64_badtype",                   field_view(4.1f),                parse_and_discard<std::int64_t>                },
        {"int64_unsigned_gtmax",            field_view(0x8000000000000000u), parse_and_discard<std::int64_t>                },
        {"int64_unsigned_absmax",           ui64_absmax,                     parse_and_discard<std::int64_t>                },

        {"uint64_null",                     field_view(),                    parse_and_discard<std::uint64_t>               },
        {"uint64_badtype",                  field_view("abc"),               parse_and_discard<std::uint64_t>               },

        {"bool_null",                       field_view(),                    parse_and_discard<bool>                        },
        {"bool_badtype",                    field_view("abc"),               parse_and_discard<bool>                        },

        {"float_null",                      field_view(),                    parse_and_discard<float>                       },
        {"float_badtype",                   field_view(42),                  parse_and_discard<float>                       },

        {"double_null",                     field_view(),                    parse_and_discard<double>                      },
        {"double_badtype",                  field_view("abc"),               parse_and_discard<double>                      },

        {"date_null",                       field_view(),                    parse_and_discard<date>                        },
        {"date_badtype",                    field_view(1.1),                 parse_and_discard<date>                        },

        {"datetime_null",                   field_view(),                    parse_and_discard<datetime>                    },
        {"datetime_badtype",                field_view("abc"),               parse_and_discard<datetime>                    },

        {"time_null",                       field_view(),                    parse_and_discard<boost::mysql::time>          },
        {"time_badtype",                    field_view(30),                  parse_and_discard<boost::mysql::time>          },

        {"string_null",                     field_view(),                    parse_and_discard<std::string>                 },
        {"string_badtype",                  field_view(makebv("abc")),       parse_and_discard<std::string>                 },

        {"blob_null",                       field_view(),                    parse_and_discard<blob>                        },
        {"blob_badtype",                    field_view("abc"),               parse_and_discard<blob>                        },
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
        {"std_optional_underlying_error",   field_view("a"),                 parse_and_discard<std::optional<std::int8_t>>  },
#endif
        {"boost_optional_underlying_error", field_view("a"),                 parse_and_discard<boost::optional<std::int8_t>>},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST(tc.parse_fn(tc.from) == client_errc::static_row_parsing_error);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
