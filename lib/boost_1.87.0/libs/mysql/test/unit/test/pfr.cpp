//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/pfr.hpp>

#if BOOST_PFR_ENABLED && defined(BOOST_MYSQL_CXX14)

#include <boost/mysql/detail/typing/pos_map.hpp>
#include <boost/mysql/detail/typing/row_traits.hpp>

#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_meta.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

using detail::get_row_name_table;
using detail::get_row_size;
using detail::is_pfr_reflectable;
using detail::is_static_row;
using detail::name_table_t;

BOOST_AUTO_TEST_SUITE(test_row_traits_pfr)

struct empty
{
};

struct s1
{
    std::int32_t i;
    float f;
    double double_field;
};

struct s2
{
    std::uint64_t s;
};

struct sbad
{
    int i;
    void* non_readable_field;
    double d;
};

union test_union
{
    int i;
    float f;
};

// is_pfr_reflectable
static_assert(is_pfr_reflectable<empty>(), "");
static_assert(is_pfr_reflectable<s1>(), "");
static_assert(is_pfr_reflectable<s2>(), "");
static_assert(!is_pfr_reflectable<const s1>(), "");
static_assert(!is_pfr_reflectable<s1&>(), "");
static_assert(!is_pfr_reflectable<const s1&>(), "");
static_assert(!is_pfr_reflectable<s1&&>(), "");
static_assert(!is_pfr_reflectable<const s1&&>(), "");
static_assert(!is_pfr_reflectable<test_union>(), "");
static_assert(!is_pfr_reflectable<s1[10]>(), "");
static_assert(!is_pfr_reflectable<int>(), "");
static_assert(!is_pfr_reflectable<const char*>(), "");
static_assert(!is_pfr_reflectable<s1*>(), "");

//
// pfr_by_name
//
#if BOOST_PFR_CORE_NAME_ENABLED

BOOST_AUTO_TEST_SUITE(pfr_by_name_)

// is_row_type
static_assert(is_static_row<pfr_by_name<empty>>, "");
static_assert(is_static_row<pfr_by_name<s1>>, "");
static_assert(is_static_row<pfr_by_name<s2>>, "");
static_assert(is_static_row<pfr_by_name<sbad>>, "");

// size
static_assert(get_row_size<pfr_by_name<empty>>() == 0u, "");
static_assert(get_row_size<pfr_by_name<s1>>() == 3u, "");
static_assert(get_row_size<pfr_by_name<s2>>() == 1u, "");

// name table
static void compare_name_tables(name_table_t lhs, name_table_t rhs)
{
    std::vector<string_view> lhsvec(lhs.begin(), lhs.end());
    std::vector<string_view> rhsvec(rhs.begin(), rhs.end());
    BOOST_TEST(lhsvec == rhsvec);
}

BOOST_AUTO_TEST_CASE(get_row_name_table_)
{
    const string_view expected_s1[] = {"i", "f", "double_field"};
    const string_view expected_s2[] = {"s"};

    compare_name_tables(get_row_name_table<pfr_by_name<empty>>(), name_table_t());
    compare_name_tables(get_row_name_table<pfr_by_name<s1>>(), expected_s1);
    compare_name_tables(get_row_name_table<pfr_by_name<s2>>(), expected_s2);
}

// meta check
BOOST_AUTO_TEST_CASE(meta_check_ok)
{
    const metadata meta[] = {
        meta_builder().type(column_type::float_).nullable(false).build(),
        meta_builder().type(column_type::double_).nullable(false).build(),
        meta_builder().type(column_type::smallint).nullable(false).build(),
    };
    const std::size_t pos_map[] = {2, 0, 1};
    diagnostics diag;
    auto err = detail::meta_check<pfr_by_name<s1>>(pos_map, meta, diag);
    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.client_message() == "");
}

BOOST_AUTO_TEST_CASE(meta_check_fail)
{
    const metadata meta[] = {
        meta_builder().type(column_type::tinyint).nullable(false).build(),
        meta_builder().type(column_type::double_).nullable(false).build(),
        meta_builder().type(column_type::double_).nullable(false).build(),
    };
    const std::size_t pos_map[] = {0, 1, 2};
    diagnostics diag;
    auto err = detail::meta_check<pfr_by_name<s1>>(pos_map, meta, diag);
    BOOST_TEST(err == client_errc::metadata_check_failed);
    BOOST_TEST(
        diag.client_message() ==
        "Incompatible types for field 'f': C++ type 'float' is not compatible with DB type 'DOUBLE'"
    );
}

BOOST_AUTO_TEST_CASE(meta_check_empty_struct)
{
    diagnostics diag;
    auto err = detail::meta_check<pfr_by_name<empty>>(
        boost::span<const std::size_t>(),
        metadata_collection_view(),
        diag
    );
    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.client_message() == "");
}

// parsing
BOOST_AUTO_TEST_CASE(parse_success)
{
    // int, float, double
    const auto fv = make_fv_arr(8.1, "abc", 42, 4.3f);
    const std::size_t pos_map[] = {2, 3, 0};
    s1 value;
    auto err = detail::parse<pfr_by_name<s1>>(pos_map, fv, value);
    BOOST_TEST(err == error_code());
    BOOST_TEST(value.i == 42);
    BOOST_TEST(value.f == 4.3f);
    BOOST_TEST(value.double_field == 8.1);
}

BOOST_AUTO_TEST_CASE(parse_error)
{
    // int, float, double
    const auto fv = make_fv_arr(8.1, "abc", nullptr, 4.3f);
    const std::size_t pos_map[] = {2, 3, 0};
    s1 value;
    auto err = detail::parse<pfr_by_name<s1>>(pos_map, fv, value);
    BOOST_TEST(err == client_errc::static_row_parsing_error);
}

BOOST_AUTO_TEST_CASE(parse_empty_struct)
{
    empty value;
    auto err = detail::parse<pfr_by_name<empty>>(
        boost::span<const std::size_t>(),
        boost::span<const field_view>(),
        value
    );
    BOOST_TEST(err == error_code());
}

BOOST_AUTO_TEST_SUITE_END()
#endif

//
// pfr_by_position
//
BOOST_AUTO_TEST_SUITE(pfr_by_position_)

// is_row_type
static_assert(is_static_row<pfr_by_position<empty>>, "");
static_assert(is_static_row<pfr_by_position<s1>>, "");
static_assert(is_static_row<pfr_by_position<s2>>, "");
static_assert(is_static_row<pfr_by_position<sbad>>, "");

// size
static_assert(get_row_size<pfr_by_position<empty>>() == 0u, "");
static_assert(get_row_size<pfr_by_position<s1>>() == 3u, "");
static_assert(get_row_size<pfr_by_position<s2>>() == 1u, "");

// name table
static_assert(get_row_name_table<pfr_by_position<empty>>().size() == 0u, "");
static_assert(get_row_name_table<pfr_by_position<s1>>().size() == 0u, "");
static_assert(get_row_name_table<pfr_by_position<s2>>().size() == 0u, "");

// meta check
BOOST_AUTO_TEST_CASE(meta_check_ok)
{
    const metadata meta[] = {
        meta_builder().type(column_type::float_).nullable(false).build(),
        meta_builder().type(column_type::double_).nullable(false).build(),
        meta_builder().type(column_type::smallint).nullable(false).build(),
    };
    const std::size_t pos_map[] = {2, 0, 1};
    diagnostics diag;
    auto err = detail::meta_check<pfr_by_position<s1>>(pos_map, meta, diag);
    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.client_message() == "");
}

BOOST_AUTO_TEST_CASE(meta_check_fail)
{
    const metadata meta[] = {
        meta_builder().type(column_type::tinyint).nullable(false).build(),
        meta_builder().type(column_type::double_).nullable(false).build(),
        meta_builder().type(column_type::double_).nullable(false).build(),
    };
    const std::size_t pos_map[] = {0, 1, 2};
    diagnostics diag;
    auto err = detail::meta_check<pfr_by_position<s1>>(pos_map, meta, diag);
    BOOST_TEST(err == client_errc::metadata_check_failed);
    BOOST_TEST(
        diag.client_message() ==
        "Incompatible types for field in position 1: C++ type 'float' is not compatible with DB type 'DOUBLE'"
    );
}

BOOST_AUTO_TEST_CASE(meta_check_empty_struct)
{
    diagnostics diag;
    auto err = detail::meta_check<pfr_by_position<empty>>(
        boost::span<const std::size_t>(),
        metadata_collection_view(),
        diag
    );
    BOOST_TEST(err == error_code());
    BOOST_TEST(diag.client_message() == "");
}

// parsing
BOOST_AUTO_TEST_CASE(parse_success)
{
    // int, float, double
    const auto fv = make_fv_arr(8.1, "abc", 42, 4.3f);
    const std::size_t pos_map[] = {2, 3, 0};
    s1 value;
    auto err = detail::parse<pfr_by_position<s1>>(pos_map, fv, value);
    BOOST_TEST(err == error_code());
    BOOST_TEST(value.i == 42);
    BOOST_TEST(value.f == 4.3f);
    BOOST_TEST(value.double_field == 8.1);
}

BOOST_AUTO_TEST_CASE(parse_error)
{
    // int, float, double
    const auto fv = make_fv_arr(8.1, "abc", nullptr, 4.3f);
    const std::size_t pos_map[] = {2, 3, 0};
    s1 value;
    auto err = detail::parse<pfr_by_position<s1>>(pos_map, fv, value);
    BOOST_TEST(err == client_errc::static_row_parsing_error);
}

BOOST_AUTO_TEST_CASE(parse_empty_struct)
{
    empty value;
    auto err = detail::parse<pfr_by_position<empty>>(
        boost::span<const std::size_t>(),
        boost::span<const field_view>(),
        value
    );
    BOOST_TEST(err == error_code());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

#endif
