//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

/**
 * These tests try to cover a range of the possible types and values MySQL support.
 * We list here all the tables that look like types_*, its contents and metadata.
 * We try reading them using the text and binary protocols, and write them using the binary
 * protocol. Table definitions must match those in db_setup.sql
 */

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/tcp.hpp>

#include <boost/mysql/detail/auxiliar/stringize.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <cstdint>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "metadata_validator.hpp"
#include "printing.hpp"
#include "safe_getenv.hpp"
#include "tcp_network_fixture.hpp"
#include "test_common.hpp"

using namespace boost::mysql::test;
using boost::mysql::blob_view;
using boost::mysql::column_type;
using boost::mysql::date;
using boost::mysql::datetime;
using boost::mysql::metadata;
using boost::mysql::results;
using boost::mysql::detail::stringize;

BOOST_AUTO_TEST_SUITE(test_database_types)

// Helpers
bool is_mariadb() { return safe_getenv("BOOST_MYSQL_TEST_DB", "mysql8") == "mariadb"; }

using flagsvec = std::vector<meta_validator::flag_getter>;

const flagsvec flags_unsigned{&metadata::is_unsigned};
const flagsvec flags_zerofill{&metadata::is_unsigned, &metadata::is_zerofill};
const flagsvec no_flags{};

struct database_types_fixture : tcp_network_fixture
{
    // Sets the time_zone to a well known value, so we can deterministically read TIMESTAMPs
    void set_time_zone()
    {
        results result;
        conn.query("SET session time_zone = '+02:00'", result);
    }

    void set_sql_mode()
    {
        results result;
        conn.query("SET session sql_mode = 'ALLOW_INVALID_DATES'", result);
    }

    database_types_fixture()
    {
        connect();
        set_time_zone();
        set_sql_mode();
    }
};

struct table
{
    std::string name;
    std::vector<meta_validator> metas;
    std::vector<boost::mysql::row> rws;

    table(std::string name) : name(std::move(name))
    {
        add_meta(
            "id",
            column_type::varchar,
            flagsvec{
                &metadata::is_primary_key,
                &metadata::is_not_null,
                &metadata::has_no_default_value,
            }
        );
    }

    void add_meta(
        std::string field,
        column_type type,
        flagsvec flags = {},
        unsigned decimals = 0,
        flagsvec ignore_flags = {}
    )
    {
        metas.emplace_back(name, std::move(field), type, std::move(flags), decimals, std::move(ignore_flags));
    }

    template <typename... Args>
    void add_row(const char* id, Args&&... args)
    {
        assert(sizeof...(Args) + 1 == metas.size());
        rws.emplace_back(makerow(id, std::forward<Args>(args)...));
    }

    void validate_rows(boost::mysql::rows_view actual_matrix) const
    {
        // The matrix size is correct
        BOOST_TEST_REQUIRE(metas.size() == actual_matrix.num_columns());

        // Build a map with the received rows, by ID
        std::unordered_map<std::string, boost::mysql::row_view> actual;
        for (const auto& row : actual_matrix)
            actual[std::string(row.at(0).as_string())] = row;

        // Verify that all expected rows are there and match
        for (const auto& expected_row : rws)
        {
            auto id = expected_row.at(0).as_string();
            BOOST_TEST_CONTEXT("row_id=" << id)
            {
                auto it = actual.find(std::string(id));
                if (it == actual.end())
                {
                    BOOST_TEST(false, "Row not found in the actual table");
                }
                else
                {
                    BOOST_TEST(
                        expected_row == it->second,
                        "\nlhs: " << expected_row << "\nrhs: " << it->second
                    );
                    actual.erase(it);
                }
            }
        }

        // Verify that there are no additional rows
        for (const auto& additional_row : actual)
        {
            BOOST_TEST_CONTEXT("row_id=" << additional_row.first)
            {
                BOOST_TEST(false, "Row was found in the table but not declared in database_types");
            }
        }
    }

    std::string select_sql() const { return stringize("SELECT * FROM ", name, " ORDER BY id"); }

    std::string insert_sql() const
    {
        std::ostringstream ss;
        ss << "INSERT INTO " << name << " VALUES (";
        for (std::size_t i = 0; i < metas.size(); ++i)
        {
            if (i == 0)
                ss << "?";
            else
                ss << ", ?";
        }
        ss << ")";
        return ss.str();
    }

    std::string delete_sql() const { return stringize("DELETE FROM ", name); }
};

// Int helpers
void int_table_columns(table& output, column_type type)
{
    output.add_meta("field_signed", type);
    output.add_meta("field_unsigned", type, flags_unsigned);
    output.add_meta("field_width", type);
    output.add_meta("field_zerofill", type, flags_zerofill);
}

table types_tinyint()
{
    table res("types_tinyint");
    int_table_columns(res, column_type::tinyint);

    res.add_row("regular", 20, 20u, 20, 20u);
    res.add_row("negative", -20, nullptr, -20, nullptr);
    res.add_row("min", -0x80, 0u, nullptr, 0u);
    res.add_row("max", 0x7f, 0xffu, nullptr, nullptr);
    return res;
}

table types_smallint()
{
    table res("types_smallint");
    int_table_columns(res, column_type::smallint);

    res.add_row("regular", 20, 20u, 20, 20u);
    res.add_row("negative", -20, nullptr, -20, nullptr);
    res.add_row("min", -0x8000, 0u, nullptr, 0u);
    res.add_row("max", 0x7fff, 0xffffu, nullptr, nullptr);
    return res;
}

table types_mediumint()
{
    table res("types_mediumint");
    int_table_columns(res, column_type::mediumint);

    res.add_row("regular", 20, 20u, 20, 20u);
    res.add_row("negative", -20, nullptr, -20, nullptr);
    res.add_row("min", -0x800000, 0u, nullptr, 0u);
    res.add_row("max", 0x7fffff, 0xffffffu, nullptr, nullptr);
    return res;
}

table types_int()
{
    table res("types_int");
    int_table_columns(res, column_type::int_);

    res.add_row("regular", 20, 20u, 20, 20u);
    res.add_row("negative", -20, nullptr, -20, nullptr);
    res.add_row("min", -0x80000000LL, 0u, nullptr, 0u);
    res.add_row("max", 0x7fffffff, 0xffffffffu, nullptr, nullptr);
    return res;
}

table types_bigint()
{
    table res("types_bigint");
    int_table_columns(res, column_type::bigint);

    res.add_row("regular", 20, 20u, 20, 20u);
    res.add_row("negative", -20, nullptr, -20, nullptr);
    res.add_row("min", -0x7fffffffffffffff - 1, 0u, nullptr, 0u);
    res.add_row("max", 0x7fffffffffffffff, 0xffffffffffffffffu, nullptr, nullptr);
    return res;
}

table types_year()
{
    table res("types_year");
    res.add_meta("field_default", column_type::year, flags_zerofill);

    res.add_row("regular", 2019u);
    res.add_row("min", 1901u);
    res.add_row("max", 2155u);
    res.add_row("zero", 0u);
    return res;
}

table types_bit()
{
    table res("types_bit");
    const char* columns[] = {
        "field_1",
        "field_8",
        "field_14",
        "field_16",
        "field_24",
        "field_25",
        "field_32",
        "field_40",
        "field_48",
        "field_56",
        "field_64"};
    for (const char* col : columns)
        res.add_meta(col, column_type::bit, flags_unsigned);

    // clang-format off
    res.add_row("min",     0u, 0u,    0u,      0u,      0u,        0u,         0u,          0u,            0u,              0u,                0u);
    res.add_row("regular", 1u, 0x9eu, 0x1e2au, 0x1234u, 0x123456u, 0x154abe0u, 0x12345678u, 0x123456789au, 0x123456789abcu, 0x123456789abcdeu, 0x1234567812345678u);
    res.add_row("max",     1u, 0xffu, 0x3fffu, 0xffffu, 0xffffffu, 0x1ffffffu, 0xffffffffu, 0xffffffffffu, 0xffffffffffffu, 0xffffffffffffffu, 0xffffffffffffffffu);
    // clang-format on
    return res;
}

table types_float()
{
    table res("types_float");
    res.add_meta("field_signed", column_type::float_, no_flags, 31);
    res.add_meta("field_unsigned", column_type::float_, flags_unsigned, 31);
    res.add_meta("field_width", column_type::float_, no_flags, 10);
    res.add_meta("field_zerofill", column_type::float_, flags_zerofill, 31);

    res.add_row("zero", 0.f, 0.f, 0.f, 0.f);
    res.add_row("int_positive", 4.f, nullptr, nullptr, nullptr);
    res.add_row("int_negative", -4.f, nullptr, nullptr, nullptr);
    res.add_row("fractional_positive", 4.2f, 4.2f, 4.2f, 4.2f);
    res.add_row("fractional_negative", -4.2f, nullptr, -4.2f, nullptr);
    res.add_row("positive_exp_positive_int", 3e20f, nullptr, nullptr, nullptr);
    res.add_row("positive_exp_negative_int", -3e20f, nullptr, nullptr, nullptr);
    res.add_row("positive_exp_positive_fractional", 3.14e20f, nullptr, nullptr, 3.14e20f);
    res.add_row("positive_exp_negative_fractional", -3.14e20f, nullptr, nullptr, nullptr);
    res.add_row("negative_exp_positive_fractional", 3.14e-20f, nullptr, nullptr, 3.14e-20f);
    return res;
}

table types_double()
{
    table res("types_double");
    res.add_meta("field_signed", column_type::double_, no_flags, 31);
    res.add_meta("field_unsigned", column_type::double_, flags_unsigned, 31);
    res.add_meta("field_width", column_type::double_, no_flags, 10);
    res.add_meta("field_zerofill", column_type::double_, flags_zerofill, 31);

    res.add_row("zero", 0., 0., 0., 0.);
    res.add_row("int_positive", 4., nullptr, nullptr, nullptr);
    res.add_row("int_negative", -4., nullptr, nullptr, nullptr);
    res.add_row("fractional_positive", 4.2, 4.2, 4.2, 4.2);
    res.add_row("fractional_negative", -4.2, nullptr, -4.2, nullptr);
    res.add_row("positive_exp_positive_int", 3e200, nullptr, nullptr, nullptr);
    res.add_row("positive_exp_negative_int", -3e200, nullptr, nullptr, nullptr);
    res.add_row("positive_exp_positive_fractional", 3.14e200, nullptr, nullptr, 3.14e200);
    res.add_row("positive_exp_negative_fractional", -3.14e200, nullptr, nullptr, nullptr);
    res.add_row("negative_exp_positive_fractional", 3.14e-200, nullptr, nullptr, 3.14e-200);
    return res;
}

table types_date()
{
    table res("types_date");
    res.add_meta("field_date", column_type::date);

    res.add_row("regular", date(2010u, 3u, 28u));
    res.add_row("leap_regular", date(1788u, 2u, 29u));
    res.add_row("leap_400", date(2000u, 2u, 29u));
    res.add_row("min", date(0u, 1u, 01u));
    res.add_row("max", date(9999u, 12u, 31u));
    res.add_row("zero", date(0u, 0u, 0u));
    res.add_row("yzero_mzero_dregular", date(0u, 0u, 20u));
    res.add_row("yzero_mregular_dzero", date(0u, 11u, 0u));
    res.add_row("yzero_invalid_date", date(0u, 11u, 31u));
    res.add_row("yregular_mzero_dzero", date(2020u, 0u, 0u));
    res.add_row("yregular_mzero_dregular", date(2020u, 0u, 20u));
    res.add_row("yregular_mregular_dzero", date(2020u, 11u, 0u));
    res.add_row("yregular_invalid_date", date(2020u, 11u, 31u));
    res.add_row("yregular_invalid_date_leapregular", date(1999u, 2u, 29u));
    res.add_row("yregular_invalid_date_leap100", date(1900u, 2u, 29u));
    return res;
}

void datetime_timestamp_common_rows(table& res)
{
    // clang-format off
    res.add_row("date",         datetime(2010, 5, 2, 0, 0, 0, 0),      datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0));
    res.add_row("date_leap4",   datetime(2004, 2, 29, 0, 0, 0, 0),     datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0));
    res.add_row("date_leap400", datetime(2000, 2, 29, 0, 0, 0, 0),     datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0));
    res.add_row("u",            nullptr,                               datetime(2010, 5, 2, 0, 0, 0, 100000),      datetime(2010, 5, 2, 0, 0, 0, 120000),      datetime(2010, 5, 2, 0, 0, 0, 123000),      datetime(2010, 5, 2, 0, 0, 0, 123400),      datetime(2010, 5, 2, 0, 0, 0, 123450),      datetime(2010, 5, 2, 0, 0, 0, 123456));
    res.add_row("s",            datetime(2010, 5, 2, 0, 0, 50, 0),     datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0)), 
    res.add_row("m",            datetime(2010, 5, 2, 0, 1, 0, 0),      datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0));
    res.add_row("hs",           datetime(2010, 5, 2, 23, 0, 50, 0),    datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0));
    res.add_row("ms",           datetime(2010, 5, 2, 0, 1, 50, 0),     datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0));
    res.add_row("hu",           nullptr,                               datetime(2010, 5, 2, 23, 0, 0, 100000),     datetime(2010, 5, 2, 23, 0, 0, 120000),     datetime(2010, 5, 2, 23, 0, 0, 123000),     datetime(2010, 5, 2, 23, 0, 0, 123400),     datetime(2010, 5, 2, 23, 0, 0, 123450),     datetime(2010, 5, 2, 23, 0, 0, 123456));
    res.add_row("mu",           nullptr,                               datetime(2010, 5, 2, 0, 1, 0, 100000),      datetime(2010, 5, 2, 0, 1, 0, 120000),      datetime(2010, 5, 2, 0, 1, 0, 123000),      datetime(2010, 5, 2, 0, 1, 0, 123400),      datetime(2010, 5, 2, 0, 1, 0, 123450),      datetime(2010, 5, 2, 0, 1, 0, 123456));
    res.add_row("hmu",          nullptr,                               datetime(2010, 5, 2, 23, 1, 0, 100000),     datetime(2010, 5, 2, 23, 1, 0, 120000),     datetime(2010, 5, 2, 23, 1, 0, 123000),     datetime(2010, 5, 2, 23, 1, 0, 123400),     datetime(2010, 5, 2, 23, 1, 0, 123450),     datetime(2010, 5, 2, 23, 1, 0, 123456));
    res.add_row("su",           nullptr,                               datetime(2010, 5, 2, 0, 0, 50, 100000),     datetime(2010, 5, 2, 0, 0, 50, 120000),     datetime(2010, 5, 2, 0, 0, 50, 123000),     datetime(2010, 5, 2, 0, 0, 50, 123400),     datetime(2010, 5, 2, 0, 0, 50, 123450),     datetime(2010, 5, 2, 0, 0, 50, 123456));
    res.add_row("hsu",          nullptr,                               datetime(2010, 5, 2, 23, 0, 50, 100000),    datetime(2010, 5, 2, 23, 0, 50, 120000),    datetime(2010, 5, 2, 23, 0, 50, 123000),    datetime(2010, 5, 2, 23, 0, 50, 123400),    datetime(2010, 5, 2, 23, 0, 50, 123450),    datetime(2010, 5, 2, 23, 0, 50, 123456));
    res.add_row("msu",          nullptr,                               datetime(2010, 5, 2, 0, 1, 50, 100000),     datetime(2010, 5, 2, 0, 1, 50, 120000),     datetime(2010, 5, 2, 0, 1, 50, 123000),     datetime(2010, 5, 2, 0, 1, 50, 123400),     datetime(2010, 5, 2, 0, 1, 50, 123450),     datetime(2010, 5, 2, 0, 1, 50, 123456));
    res.add_row("h",            datetime(2010, 5, 2, 23, 0, 0, 0),     datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0));
    res.add_row("hm",           datetime(2010, 5, 2, 23, 1, 0, 0),     datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0));
    res.add_row("hms",          datetime(2010, 5, 2, 23, 1, 50, 0),    datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0));
    res.add_row("hmsu",         nullptr,                               datetime(2010, 5, 2, 23, 1, 50, 100000),    datetime(2010, 5, 2, 23, 1, 50, 120000),    datetime(2010, 5, 2, 23, 1, 50, 123000),    datetime(2010, 5, 2, 23, 1, 50, 123400),    datetime(2010, 5, 2, 23, 1, 50, 123450),    datetime(2010, 5, 2, 23, 1, 50, 123456));
    // clang-format on
}

table types_datetime()
{
    table res("types_datetime");
    res.add_meta("field_0", column_type::datetime, no_flags, 0, flags_unsigned);
    res.add_meta("field_1", column_type::datetime, no_flags, 1, flags_unsigned);
    res.add_meta("field_2", column_type::datetime, no_flags, 2, flags_unsigned);
    res.add_meta("field_3", column_type::datetime, no_flags, 3, flags_unsigned);
    res.add_meta("field_4", column_type::datetime, no_flags, 4, flags_unsigned);
    res.add_meta("field_5", column_type::datetime, no_flags, 5, flags_unsigned);
    res.add_meta("field_6", column_type::datetime, no_flags, 6, flags_unsigned);

    datetime_timestamp_common_rows(res);

    // clang-format off
    res.add_row("min",          datetime(0, 1, 1, 0, 0, 0, 0),         datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0));
    res.add_row("max",          datetime(9999, 12, 31, 23, 59, 59, 0), datetime(9999, 12, 31, 23, 59, 59, 900000), datetime(9999, 12, 31, 23, 59, 59, 990000), datetime(9999, 12, 31, 23, 59, 59, 999000), datetime(9999, 12, 31, 23, 59, 59, 999900), datetime(9999, 12, 31, 23, 59, 59, 999990), datetime(9999, 12, 31, 23, 59, 59, 999999));
    
    res.add_row("date_zero",                               datetime(   0,  0,  0,  0,  0,  0, 0), datetime(   0,  0,  0,  0,  0,  0, 0), datetime(   0,  0,  0,  0,  0,  0,  0), datetime(   0,  0,  0,  0,  0,  0,   0), datetime(   0,  0,  0,  0,  0,  0,    0), datetime(   0,  0,  0,  0,  0,  0,     0), datetime(   0,  0,  0,  0,  0,  0,      0));
    res.add_row("date_yzero_mzero_dregular",               datetime(   0,  0, 10,  0,  0,  0, 0), datetime(   0,  0, 10,  0,  0,  0, 0), datetime(   0,  0, 10,  0,  0,  0,  0), datetime(   0,  0, 10,  0,  0,  0,   0), datetime(   0,  0, 10,  0,  0,  0,    0), datetime(   0,  0, 10,  0,  0,  0,     0), datetime(   0,  0, 10,  0,  0,  0,      0));
    res.add_row("date_yzero_mregular_dzero",               datetime(   0, 10,  0,  0,  0,  0, 0), datetime(   0, 10,  0,  0,  0,  0, 0), datetime(   0, 10,  0,  0,  0,  0,  0), datetime(   0, 10,  0,  0,  0,  0,   0), datetime(   0, 10,  0,  0,  0,  0,    0), datetime(   0, 10,  0,  0,  0,  0,     0), datetime(   0, 10,  0,  0,  0,  0,      0));
    res.add_row("date_yzero_invalid_date",                 datetime(   0, 11, 31,  0,  0,  0, 0), datetime(   0, 11, 31,  0,  0,  0, 0), datetime(   0, 11, 31,  0,  0,  0,  0), datetime(   0, 11, 31,  0,  0,  0,   0), datetime(   0, 11, 31,  0,  0,  0,    0), datetime(   0, 11, 31,  0,  0,  0,     0), datetime(   0, 11, 31,  0,  0,  0,      0));
    res.add_row("date_yregular_mzero_dzero",               datetime(2020,  0,  0,  0,  0,  0, 0), datetime(2020,  0,  0,  0,  0,  0, 0), datetime(2020,  0,  0,  0,  0,  0,  0), datetime(2020,  0,  0,  0,  0,  0,   0), datetime(2020,  0,  0,  0,  0,  0,    0), datetime(2020,  0,  0,  0,  0,  0,     0), datetime(2020,  0,  0,  0,  0,  0,      0));
    res.add_row("date_yregular_mzero_dregular",            datetime(2020,  0, 10,  0,  0,  0, 0), datetime(2020,  0, 10,  0,  0,  0, 0), datetime(2020,  0, 10,  0,  0,  0,  0), datetime(2020,  0, 10,  0,  0,  0,   0), datetime(2020,  0, 10,  0,  0,  0,    0), datetime(2020,  0, 10,  0,  0,  0,     0), datetime(2020,  0, 10,  0,  0,  0,      0));
    res.add_row("date_yregular_mregular_dzero",            datetime(2020, 10,  0,  0,  0,  0, 0), datetime(2020, 10,  0,  0,  0,  0, 0), datetime(2020, 10,  0,  0,  0,  0,  0), datetime(2020, 10,  0,  0,  0,  0,   0), datetime(2020, 10,  0,  0,  0,  0,    0), datetime(2020, 10,  0,  0,  0,  0,     0), datetime(2020, 10,  0,  0,  0,  0,      0));
    res.add_row("date_yregular_invalid_date",              datetime(2020, 11, 31,  0,  0,  0, 0), datetime(2020, 11, 31,  0,  0,  0, 0), datetime(2020, 11, 31,  0,  0,  0,  0), datetime(2020, 11, 31,  0,  0,  0,   0), datetime(2020, 11, 31,  0,  0,  0,    0), datetime(2020, 11, 31,  0,  0,  0,     0), datetime(2020, 11, 31,  0,  0,  0,      0));
    res.add_row("date_yregular_invalid_date_leapregular",  datetime(1999,  2, 29,  0,  0,  0, 0), datetime(1999,  2, 29,  0,  0,  0, 0), datetime(1999,  2, 29,  0,  0,  0,  0), datetime(1999,  2, 29,  0,  0,  0,   0), datetime(1999,  2, 29,  0,  0,  0,    0), datetime(1999,  2, 29,  0,  0,  0,     0), datetime(1999,  2, 29,  0,  0,  0,      0));
    res.add_row("date_yregular_invalid_date_leap100",      datetime(1900,  2, 29,  0,  0,  0, 0), datetime(1900,  2, 29,  0,  0,  0, 0), datetime(1900,  2, 29,  0,  0,  0,  0), datetime(1900,  2, 29,  0,  0,  0,   0), datetime(1900,  2, 29,  0,  0,  0,    0), datetime(1900,  2, 29,  0,  0,  0,     0), datetime(1900,  2, 29,  0,  0,  0,      0)),
    
    res.add_row("hms_zero",                              datetime(   0,  0,  0, 10, 20, 30, 0), datetime(   0,  0,  0, 10, 20, 30, 0), datetime(   0,  0,  0, 10, 20, 30,  0), datetime(   0,  0,  0, 10, 20, 30,   0), datetime(   0,  0,  0, 10, 20, 30,    0), datetime(   0,  0,  0, 10, 20, 30,     0), datetime(   0,  0,  0, 10, 20, 30,      0));
    res.add_row("hms_yzero_mzero_dregular",              datetime(   0,  0, 10, 10, 20, 30, 0), datetime(   0,  0, 10, 10, 20, 30, 0), datetime(   0,  0, 10, 10, 20, 30,  0), datetime(   0,  0, 10, 10, 20, 30,   0), datetime(   0,  0, 10, 10, 20, 30,    0), datetime(   0,  0, 10, 10, 20, 30,     0), datetime(   0,  0, 10, 10, 20, 30,      0));
    res.add_row("hms_yzero_mregular_dzero",              datetime(   0, 10,  0, 10, 20, 30, 0), datetime(   0, 10,  0, 10, 20, 30, 0), datetime(   0, 10,  0, 10, 20, 30,  0), datetime(   0, 10,  0, 10, 20, 30,   0), datetime(   0, 10,  0, 10, 20, 30,    0), datetime(   0, 10,  0, 10, 20, 30,     0), datetime(   0, 10,  0, 10, 20, 30,      0));
    res.add_row("hms_yzero_invalid_date",                datetime(   0, 11, 31, 10, 20, 30, 0), datetime(   0, 11, 31, 10, 20, 30, 0), datetime(   0, 11, 31, 10, 20, 30,  0), datetime(   0, 11, 31, 10, 20, 30,   0), datetime(   0, 11, 31, 10, 20, 30,    0), datetime(   0, 11, 31, 10, 20, 30,     0), datetime(   0, 11, 31, 10, 20, 30,      0));
    res.add_row("hms_yregular_mzero_dzero",              datetime(2020,  0,  0, 10, 20, 30, 0), datetime(2020,  0,  0, 10, 20, 30, 0), datetime(2020,  0,  0, 10, 20, 30,  0), datetime(2020,  0,  0, 10, 20, 30,   0), datetime(2020,  0,  0, 10, 20, 30,    0), datetime(2020,  0,  0, 10, 20, 30,     0), datetime(2020,  0,  0, 10, 20, 30,      0));
    res.add_row("hms_yregular_mzero_dregular",           datetime(2020,  0, 10, 10, 20, 30, 0), datetime(2020,  0, 10, 10, 20, 30, 0), datetime(2020,  0, 10, 10, 20, 30,  0), datetime(2020,  0, 10, 10, 20, 30,   0), datetime(2020,  0, 10, 10, 20, 30,    0), datetime(2020,  0, 10, 10, 20, 30,     0), datetime(2020,  0, 10, 10, 20, 30,      0));
    res.add_row("hms_yregular_mregular_dzero",           datetime(2020, 10,  0, 10, 20, 30, 0), datetime(2020, 10,  0, 10, 20, 30, 0), datetime(2020, 10,  0, 10, 20, 30,  0), datetime(2020, 10,  0, 10, 20, 30,   0), datetime(2020, 10,  0, 10, 20, 30,    0), datetime(2020, 10,  0, 10, 20, 30,     0), datetime(2020, 10,  0, 10, 20, 30,      0));
    res.add_row("hms_yregular_invalid_date",             datetime(2020, 11, 31, 10, 20, 30, 0), datetime(2020, 11, 31, 10, 20, 30, 0), datetime(2020, 11, 31, 10, 20, 30,  0), datetime(2020, 11, 31, 10, 20, 30,   0), datetime(2020, 11, 31, 10, 20, 30,    0), datetime(2020, 11, 31, 10, 20, 30,     0), datetime(2020, 11, 31, 10, 20, 30,      0));
    res.add_row("hms_yregular_invalid_date_leapregular", datetime(1999,  2, 29, 10, 20, 30, 0), datetime(1999,  2, 29, 10, 20, 30, 0), datetime(1999,  2, 29, 10, 20, 30,  0), datetime(1999,  2, 29, 10, 20, 30,   0), datetime(1999,  2, 29, 10, 20, 30,    0), datetime(1999,  2, 29, 10, 20, 30,     0), datetime(1999,  2, 29, 10, 20, 30,      0));
    res.add_row("hms_yregular_invalid_date_leap100",     datetime(1900,  2, 29, 10, 20, 30, 0), datetime(1900,  2, 29, 10, 20, 30, 0), datetime(1900,  2, 29, 10, 20, 30,  0), datetime(1900,  2, 29, 10, 20, 30,   0), datetime(1900,  2, 29, 10, 20, 30,    0), datetime(1900,  2, 29, 10, 20, 30,     0), datetime(1900,  2, 29, 10, 20, 30,      0)),
    
    
    res.add_row("hmsu_zero",                              datetime(   0,  0,  0, 10, 20, 30, 0), datetime(   0,  0,  0, 10, 20, 30, 900000), datetime(   0,  0,  0, 10, 20, 30, 990000), datetime(   0,  0,  0, 10, 20, 30, 999000), datetime(   0,  0,  0, 10, 20, 30, 999900), datetime(   0,  0,  0, 10, 20, 30, 999990), datetime(   0,  0,  0, 10, 20, 30, 999999));
    res.add_row("hmsu_yzero_mzero_dregular",              datetime(   0,  0, 10, 10, 20, 30, 0), datetime(   0,  0, 10, 10, 20, 30, 900000), datetime(   0,  0, 10, 10, 20, 30, 990000), datetime(   0,  0, 10, 10, 20, 30, 999000), datetime(   0,  0, 10, 10, 20, 30, 999900), datetime(   0,  0, 10, 10, 20, 30, 999990), datetime(   0,  0, 10, 10, 20, 30, 999999));
    res.add_row("hmsu_yzero_mregular_dzero",              datetime(   0, 10,  0, 10, 20, 30, 0), datetime(   0, 10,  0, 10, 20, 30, 900000), datetime(   0, 10,  0, 10, 20, 30, 990000), datetime(   0, 10,  0, 10, 20, 30, 999000), datetime(   0, 10,  0, 10, 20, 30, 999900), datetime(   0, 10,  0, 10, 20, 30, 999990), datetime(   0, 10,  0, 10, 20, 30, 999999));
    res.add_row("hmsu_yzero_invalid_date",                datetime(   0, 11, 31, 10, 20, 30, 0), datetime(   0, 11, 31, 10, 20, 30, 900000), datetime(   0, 11, 31, 10, 20, 30, 990000), datetime(   0, 11, 31, 10, 20, 30, 999000), datetime(   0, 11, 31, 10, 20, 30, 999900), datetime(   0, 11, 31, 10, 20, 30, 999990), datetime(   0, 11, 31, 10, 20, 30, 999999));
    res.add_row("hmsu_yregular_mzero_dzero",              datetime(2020,  0,  0, 10, 20, 30, 0), datetime(2020,  0,  0, 10, 20, 30, 900000), datetime(2020,  0,  0, 10, 20, 30, 990000), datetime(2020,  0,  0, 10, 20, 30, 999000), datetime(2020,  0,  0, 10, 20, 30, 999900), datetime(2020,  0,  0, 10, 20, 30, 999990), datetime(2020,  0,  0, 10, 20, 30, 999999));
    res.add_row("hmsu_yregular_mzero_dregular",           datetime(2020,  0, 10, 10, 20, 30, 0), datetime(2020,  0, 10, 10, 20, 30, 900000), datetime(2020,  0, 10, 10, 20, 30, 990000), datetime(2020,  0, 10, 10, 20, 30, 999000), datetime(2020,  0, 10, 10, 20, 30, 999900), datetime(2020,  0, 10, 10, 20, 30, 999990), datetime(2020,  0, 10, 10, 20, 30, 999999));
    res.add_row("hmsu_yregular_mregular_dzero",           datetime(2020, 10,  0, 10, 20, 30, 0), datetime(2020, 10,  0, 10, 20, 30, 900000), datetime(2020, 10,  0, 10, 20, 30, 990000), datetime(2020, 10,  0, 10, 20, 30, 999000), datetime(2020, 10,  0, 10, 20, 30, 999900), datetime(2020, 10,  0, 10, 20, 30, 999990), datetime(2020, 10,  0, 10, 20, 30, 999999));
    res.add_row("hmsu_yregular_invalid_date",             datetime(2020, 11, 31, 10, 20, 30, 0), datetime(2020, 11, 31, 10, 20, 30, 900000), datetime(2020, 11, 31, 10, 20, 30, 990000), datetime(2020, 11, 31, 10, 20, 30, 999000), datetime(2020, 11, 31, 10, 20, 30, 999900), datetime(2020, 11, 31, 10, 20, 30, 999990), datetime(2020, 11, 31, 10, 20, 30, 999999));
    res.add_row("hmsu_yregular_invalid_date_leapregular", datetime(1999,  2, 29, 10, 20, 30, 0), datetime(1999,  2, 29, 10, 20, 30, 900000), datetime(1999,  2, 29, 10, 20, 30, 990000), datetime(1999,  2, 29, 10, 20, 30, 999000), datetime(1999,  2, 29, 10, 20, 30, 999900), datetime(1999,  2, 29, 10, 20, 30, 999990), datetime(1999,  2, 29, 10, 20, 30, 999999));
    res.add_row("hmsu_yregular_invalid_date_leap100",     datetime(1900,  2, 29, 10, 20, 30, 0), datetime(1900,  2, 29, 10, 20, 30, 900000), datetime(1900,  2, 29, 10, 20, 30, 990000), datetime(1900,  2, 29, 10, 20, 30, 999000), datetime(1900,  2, 29, 10, 20, 30, 999900), datetime(1900,  2, 29, 10, 20, 30, 999990), datetime(1900,  2, 29, 10, 20, 30, 999999));
    // clang-format on
    return res;
}

table types_timestamp()
{
    table res("types_timestamp");
    res.add_meta("field_0", column_type::timestamp, no_flags, 0, flags_unsigned);
    res.add_meta("field_1", column_type::timestamp, no_flags, 1, flags_unsigned);
    res.add_meta("field_2", column_type::timestamp, no_flags, 2, flags_unsigned);
    res.add_meta("field_3", column_type::timestamp, no_flags, 3, flags_unsigned);
    res.add_meta("field_4", column_type::timestamp, no_flags, 4, flags_unsigned);
    res.add_meta("field_5", column_type::timestamp, no_flags, 5, flags_unsigned);
    res.add_meta("field_6", column_type::timestamp, no_flags, 6, flags_unsigned);

    datetime_timestamp_common_rows(res);

    // clang-format off
    res.add_row("zero", datetime(),                            datetime(),                                 datetime(),                                 datetime(),                                 datetime(),                                 datetime(),                                 datetime());
    res.add_row("min",  datetime(1970,  1,  1,  2,  0,  1, 0), datetime(1970,  1,  1,  2,  0,  1, 0),      datetime(1970,  1,  1,  2,  0,  1,  0),     datetime(1970,  1,  1,  2,  0,  1,   0),    datetime(1970,  1,  1,  2,  0,  1,    0),   datetime(1970,  1,  1,  2,  0,  1,     0),  datetime(1970,  1,  1,  2,  0,  1,      0));
    res.add_row("max",  datetime(2038,  1, 19,  5, 14,  7, 0), datetime(2038,  1, 19,  5, 14,  7, 900000), datetime(2038,  1, 19,  5, 14,  7, 990000), datetime(2038,  1, 19,  5, 14,  7, 999000), datetime(2038,  1, 19,  5, 14,  7, 999900), datetime(2038,  1, 19,  5, 14,  7, 999990), datetime(2038,  1, 19,  5, 14,  7, 999999));
    // clang-format on
    return res;
}

table types_time()
{
    table res("types_time");
    res.add_meta("field_0", column_type::time, no_flags, 0, flags_unsigned);
    res.add_meta("field_1", column_type::time, no_flags, 1, flags_unsigned);
    res.add_meta("field_2", column_type::time, no_flags, 2, flags_unsigned);
    res.add_meta("field_3", column_type::time, no_flags, 3, flags_unsigned);
    res.add_meta("field_4", column_type::time, no_flags, 4, flags_unsigned);
    res.add_meta("field_5", column_type::time, no_flags, 5, flags_unsigned);
    res.add_meta("field_6", column_type::time, no_flags, 6, flags_unsigned);

    // clang-format off
    res.add_row("zero",            maket( 0,  0,  0),  maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0));
    res.add_row("d",               maket(48,  0,  0),  maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0));
    res.add_row("negative_d",     -maket(48,  0,  0), -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0));
    res.add_row("h",               maket(23,  0,  0),  maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0));
    res.add_row("negative_h",     -maket(23,  0,  0), -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0));
    res.add_row("dh",              maket(71,  0,  0),  maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0));
    res.add_row("negative_dh",    -maket(71,  0,  0), -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0));
    res.add_row("m",               maket( 0,  1,  0),  maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0));
    res.add_row("negative_m",     -maket( 0,  1,  0), -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0));
    res.add_row("dm",              maket(48,  1,  0),  maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0));
    res.add_row("negative_dm",    -maket(48,  1,  0), -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0));
    res.add_row("hm",              maket(23,  1,  0),  maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0));
    res.add_row("negative_hm",    -maket(23,  1,  0), -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0));
    res.add_row("dhm",             maket(71,  1,  0),  maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0));
    res.add_row("negative_dhm",   -maket(71,  1,  0), -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0));
    res.add_row("s",               maket( 0,  0, 50),  maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50));
    res.add_row("negative_s",     -maket( 0,  0, 50), -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50));
    res.add_row("ds",              maket(48,  0, 50),  maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50));
    res.add_row("negative_ds",    -maket(48,  0, 50), -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50));
    res.add_row("hs",              maket(23,  0, 50),  maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50));
    res.add_row("negative_hs",    -maket(23,  0, 50), -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50));
    res.add_row("dhs",             maket(71,  0, 50),  maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50));
    res.add_row("negative_dhs",   -maket(71,  0, 50), -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50));
    res.add_row("ms",              maket( 0,  1, 50),  maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50));
    res.add_row("negative_ms",    -maket( 0,  1, 50), -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50));
    res.add_row("dms",             maket(48,  1, 50),  maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50));
    res.add_row("negative_dms",   -maket(48,  1, 50), -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50));
    res.add_row("hms",             maket(23,  1, 50),  maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50));
    res.add_row("negative_hms",   -maket(23,  1, 50), -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50));
    res.add_row("dhms",            maket(71,  1, 50),  maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50));
    res.add_row("negative_dhms",  -maket(71,  1, 50), -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50));
    res.add_row("u",              nullptr,             maket( 0,  0,  0, 100000),  maket( 0,  0,  0, 120000),  maket( 0,  0,  0, 123000),  maket( 0,  0,  0, 123400),  maket( 0,  0,  0, 123450),  maket( 0,  0,  0,  123456));
    res.add_row("negative_u",     nullptr,            -maket( 0,  0,  0, 100000), -maket( 0,  0,  0, 120000), -maket( 0,  0,  0, 123000), -maket( 0,  0,  0, 123400), -maket( 0,  0,  0, 123450), -maket( 0,  0,  0,  123456));
    res.add_row("du",             nullptr,             maket(48,  0,  0, 100000),  maket(48,  0,  0, 120000),  maket(48,  0,  0, 123000),  maket(48,  0,  0, 123400),  maket(48,  0,  0, 123450),  maket(48,  0,  0,  123456));
    res.add_row("negative_du",    nullptr,            -maket(48,  0,  0, 100000), -maket(48,  0,  0, 120000), -maket(48,  0,  0, 123000), -maket(48,  0,  0, 123400), -maket(48,  0,  0, 123450), -maket(48,  0,  0,  123456));
    res.add_row("hu",             nullptr,             maket(23,  0,  0, 100000),  maket(23,  0,  0, 120000),  maket(23,  0,  0, 123000),  maket(23,  0,  0, 123400),  maket(23,  0,  0, 123450),  maket(23,  0,  0,  123456));
    res.add_row("negative_hu",    nullptr,            -maket(23,  0,  0, 100000), -maket(23,  0,  0, 120000), -maket(23,  0,  0, 123000), -maket(23,  0,  0, 123400), -maket(23,  0,  0, 123450), -maket(23,  0,  0,  123456));
    res.add_row("dhu",            nullptr,             maket(71,  0,  0, 100000),  maket(71,  0,  0, 120000),  maket(71,  0,  0, 123000),  maket(71,  0,  0, 123400),  maket(71,  0,  0, 123450),  maket(71,  0,  0,  123456));
    res.add_row("negative_dhu",   nullptr,            -maket(71,  0,  0, 100000), -maket(71,  0,  0, 120000), -maket(71,  0,  0, 123000), -maket(71,  0,  0, 123400), -maket(71,  0,  0, 123450), -maket(71,  0,  0,  123456));
    res.add_row("mu",             nullptr,             maket( 0,  1,  0, 100000),  maket( 0,  1,  0, 120000),  maket( 0,  1,  0, 123000),  maket( 0,  1,  0, 123400),  maket( 0,  1,  0, 123450),  maket( 0,  1,  0,  123456));
    res.add_row("negative_mu",    nullptr,            -maket( 0,  1,  0, 100000), -maket( 0,  1,  0, 120000), -maket( 0,  1,  0, 123000), -maket( 0,  1,  0, 123400), -maket( 0,  1,  0, 123450), -maket( 0,  1,  0,  123456));
    res.add_row("dmu",            nullptr,             maket(48,  1,  0, 100000),  maket(48,  1,  0, 120000),  maket(48,  1,  0, 123000),  maket(48,  1,  0, 123400),  maket(48,  1,  0, 123450),  maket(48,  1,  0,  123456));
    res.add_row("negative_dmu",   nullptr,            -maket(48,  1,  0, 100000), -maket(48,  1,  0, 120000), -maket(48,  1,  0, 123000), -maket(48,  1,  0, 123400), -maket(48,  1,  0, 123450), -maket(48,  1,  0,  123456));
    res.add_row("hmu",            nullptr,             maket(23,  1,  0, 100000),  maket(23,  1,  0, 120000),  maket(23,  1,  0, 123000),  maket(23,  1,  0, 123400),  maket(23,  1,  0, 123450),  maket(23,  1,  0,  123456));
    res.add_row("negative_hmu",   nullptr,            -maket(23,  1,  0, 100000), -maket(23,  1,  0, 120000), -maket(23,  1,  0, 123000), -maket(23,  1,  0, 123400), -maket(23,  1,  0, 123450), -maket(23,  1,  0,  123456));
    res.add_row("dhmu",           nullptr,             maket(71,  1,  0, 100000),  maket(71,  1,  0, 120000),  maket(71,  1,  0, 123000),  maket(71,  1,  0, 123400),  maket(71,  1,  0, 123450),  maket(71,  1,  0,  123456));
    res.add_row("negative_dhmu",  nullptr,            -maket(71,  1,  0, 100000), -maket(71,  1,  0, 120000), -maket(71,  1,  0, 123000), -maket(71,  1,  0, 123400), -maket(71,  1,  0, 123450), -maket(71,  1,  0,  123456));
    res.add_row("su",             nullptr,             maket( 0,  0, 50, 100000),  maket( 0,  0, 50, 120000),  maket( 0,  0, 50, 123000),  maket( 0,  0, 50, 123400),  maket( 0,  0, 50, 123450),  maket( 0,  0, 50,  123456));
    res.add_row("negative_su",    nullptr,            -maket( 0,  0, 50, 100000), -maket( 0,  0, 50, 120000), -maket( 0,  0, 50, 123000), -maket( 0,  0, 50, 123400), -maket( 0,  0, 50, 123450), -maket( 0,  0, 50,  123456));
    res.add_row("dsu",            nullptr,             maket(48,  0, 50, 100000),  maket(48,  0, 50, 120000),  maket(48,  0, 50, 123000),  maket(48,  0, 50, 123400),  maket(48,  0, 50, 123450),  maket(48,  0, 50,  123456));
    res.add_row("negative_dsu",   nullptr,            -maket(48,  0, 50, 100000), -maket(48,  0, 50, 120000), -maket(48,  0, 50, 123000), -maket(48,  0, 50, 123400), -maket(48,  0, 50, 123450), -maket(48,  0, 50,  123456));
    res.add_row("hsu",            nullptr,             maket(23,  0, 50, 100000),  maket(23,  0, 50, 120000),  maket(23,  0, 50, 123000),  maket(23,  0, 50, 123400),  maket(23,  0, 50, 123450),  maket(23,  0, 50,  123456));
    res.add_row("negative_hsu",   nullptr,            -maket(23,  0, 50, 100000), -maket(23,  0, 50, 120000), -maket(23,  0, 50, 123000), -maket(23,  0, 50, 123400), -maket(23,  0, 50, 123450), -maket(23,  0, 50,  123456));
    res.add_row("dhsu",           nullptr,             maket(71,  0, 50, 100000),  maket(71,  0, 50, 120000),  maket(71,  0, 50, 123000),  maket(71,  0, 50, 123400),  maket(71,  0, 50, 123450),  maket(71,  0, 50,  123456));
    res.add_row("negative_dhsu",  nullptr,            -maket(71,  0, 50, 100000), -maket(71,  0, 50, 120000), -maket(71,  0, 50, 123000), -maket(71,  0, 50, 123400), -maket(71,  0, 50, 123450), -maket(71,  0, 50,  123456));
    res.add_row("msu",            nullptr,             maket( 0,  1, 50, 100000),  maket( 0,  1, 50, 120000),  maket( 0,  1, 50, 123000),  maket( 0,  1, 50, 123400),  maket( 0,  1, 50, 123450),  maket( 0,  1, 50,  123456));
    res.add_row("negative_msu",   nullptr,            -maket( 0,  1, 50, 100000), -maket( 0,  1, 50, 120000), -maket( 0,  1, 50, 123000), -maket( 0,  1, 50, 123400), -maket( 0,  1, 50, 123450), -maket( 0,  1, 50,  123456));
    res.add_row("dmsu",           nullptr,             maket(48,  1, 50, 100000),  maket(48,  1, 50, 120000),  maket(48,  1, 50, 123000),  maket(48,  1, 50, 123400),  maket(48,  1, 50, 123450),  maket(48,  1, 50,  123456));
    res.add_row("negative_dmsu",  nullptr,            -maket(48,  1, 50, 100000), -maket(48,  1, 50, 120000), -maket(48,  1, 50, 123000), -maket(48,  1, 50, 123400), -maket(48,  1, 50, 123450), -maket(48,  1, 50,  123456));
    res.add_row("hmsu",           nullptr,             maket(23,  1, 50, 100000),  maket(23,  1, 50, 120000),  maket(23,  1, 50, 123000),  maket(23,  1, 50, 123400),  maket(23,  1, 50, 123450),  maket(23,  1, 50,  123456));
    res.add_row("negative_hmsu",  nullptr,            -maket(23,  1, 50, 100000), -maket(23,  1, 50, 120000), -maket(23,  1, 50, 123000), -maket(23,  1, 50, 123400), -maket(23,  1, 50, 123450), -maket(23,  1, 50,  123456));
    res.add_row("dhmsu",          nullptr,             maket(71,  1, 50, 100000),  maket(71,  1, 50, 120000),  maket(71,  1, 50, 123000),  maket(71,  1, 50, 123400),  maket(71,  1, 50, 123450),  maket(71,  1, 50,  123456));
    res.add_row("negative_dhmsu", nullptr,            -maket(71,  1, 50, 100000), -maket(71,  1, 50, 120000), -maket(71,  1, 50, 123000), -maket(71,  1, 50, 123400), -maket(71,  1, 50, 123450), -maket(71,  1, 50,  123456));
    res.add_row("min",            -maket(838, 59, 59),-maket(838,59, 58, 900000), -maket(838, 59, 58, 990000),-maket(838, 59, 58, 999000),-maket(838, 59, 58, 999900),-maket(838,59, 58, 999990), -maket(838,59, 58,  999999));
    res.add_row("max",             maket(838, 59, 59), maket(838,59, 58, 900000),  maket(838, 59, 58, 990000), maket(838, 59, 58, 999000), maket(838, 59, 58, 999900), maket(838,59, 58, 999990),  maket(838,59, 58,  999999));
    // clang-format on
    return res;
}

table types_string()
{
    table res("types_string");
    res.add_meta("field_char", column_type::char_);
    res.add_meta("field_varchar", column_type::varchar);
    res.add_meta("field_tinytext", column_type::text);
    res.add_meta("field_text", column_type::text);
    res.add_meta("field_mediumtext", column_type::text);
    res.add_meta("field_longtext", column_type::text);
    res.add_meta("field_text_bincol", column_type::text);
    res.add_meta("field_enum", column_type::enum_);
    res.add_meta("field_set", column_type::set);

    // clang-format off
    res.add_row("regular", "test_char", "test_varchar", "test_tinytext", "test_text", "test_mediumtext", "test_longtext", "test_bincol", "red",    "red,green");
    res.add_row("utf8",    "\xc3\xb1",  "\xc3\x91",     "\xc3\xa1",      "\xc3\xa9",  "\xc3\xad",        "\xc3\xb3",      "\xc3\xba",    nullptr,  nullptr);
    res.add_row("empty",   "",          "",             "",              "",          "",                "",              "",            nullptr,  "");
    // clang-format on
    return res;
}

// MariaDB doesn't have a dedicated column type, so there is a difference in metadata.
// Values should be the same, though.
table types_json()
{
    table res("types_json");
    res.add_meta("field_json", is_mariadb() ? column_type::text : column_type::json);

    res.add_row("regular", R"([null, 42, false, "abc", {"key": "value"}])");
    res.add_row("unicode_escape", R"(["\\u0000value\\u0000"])");
    res.add_row("utf8", "[\"adi\xc3\xb3os\"]");
    res.add_row("empty", "{}");
    return res;
}

table types_binary()
{
    table res("types_binary");
    res.add_meta("field_binary", column_type::binary);
    res.add_meta("field_varbinary", column_type::varbinary);
    res.add_meta("field_tinyblob", column_type::blob);
    res.add_meta("field_blob", column_type::blob);
    res.add_meta("field_mediumblob", column_type::blob);
    res.add_meta("field_longblob", column_type::blob);

    // clang-format off
    res.add_row("regular",  makebv("\0_binary\0\0"),          makebv("\0_varbinary"),  makebv("\0_tinyblob"),  makebv("\0_blob"),  makebv("\0_mediumblob"), makebv("\0_longblob"));
    res.add_row("nonascii", makebv("\0\xff\0\0\0\0\0\0\0\0"), makebv("\1\xfe"),        makebv("\2\xfd"),       makebv("\3\xfc"),   makebv("\4\xfb"),        makebv("\5\xfa"));
    res.add_row("empty",    makebv("\0\0\0\0\0\0\0\0\0\0"),   blob_view(),             blob_view(),            blob_view(),        blob_view(),             blob_view());
    // clang-format on
    return res;
}

// These types don't have a better representation, and we represent
// them as strings or binary
table types_not_implemented()
{
    table res("types_not_implemented");
    res.add_meta("field_decimal", column_type::decimal);
    res.add_meta("field_geometry", column_type::geometry);

    static constexpr std::uint8_t geometry_value[] = {0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00,
                                                      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f, 0x00,
                                                      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40};

    res.add_row("regular", "300", blob_view(geometry_value));
    return res;
}

// Tests for certain metadata flags
table types_flags()
{
    table res("types_flags");
    res.add_meta(
        "field_timestamp",
        column_type::timestamp,
        flagsvec{&metadata::is_set_to_now_on_update},
        0,
        flags_unsigned
    );
    res.add_meta(
        "field_primary_key",
        column_type::int_,
        flagsvec{&metadata::is_primary_key, &metadata::is_not_null, &metadata::is_auto_increment}
    );
    res.add_meta("field_not_null", column_type::char_, flagsvec{&metadata::is_not_null});
    res.add_meta("field_unique", column_type::int_, flagsvec{&metadata::is_unique_key});
    res.add_meta("field_indexed", column_type::int_, flagsvec{&metadata::is_multiple_key});

    res.add_row("default", nullptr, 50, "char", 21, 42);
    return res;
}

std::vector<table> make_all_tables()
{
    std::vector<table> res;
    res.push_back(types_tinyint());
    res.push_back(types_smallint());
    res.push_back(types_mediumint());
    res.push_back(types_int());
    res.push_back(types_bigint());
    res.push_back(types_year());
    res.push_back(types_bit());
    res.push_back(types_float());
    res.push_back(types_double());
    res.push_back(types_date());
    res.push_back(types_datetime());
    res.push_back(types_timestamp());
    res.push_back(types_time());
    res.push_back(types_string());
    res.push_back(types_binary());
    res.push_back(types_not_implemented());
    res.push_back(types_flags());
    return res;
}

const std::vector<table>& all_tables()
{
    static std::vector<table> res = make_all_tables();
    return res;
}

BOOST_FIXTURE_TEST_CASE(query_read, database_types_fixture)
{
    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table.name)
        {
            // Execute the query
            results result;
            conn.query(table.select_sql(), result);

            // Validate the received contents
            validate_meta(result.meta(), table.metas);
            table.validate_rows(result.rows());
        }
    }
}

BOOST_FIXTURE_TEST_CASE(statement_read, database_types_fixture)
{
    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table.name)
        {
            // Prepare the statement
            auto stmt = conn.prepare_statement(table.select_sql());

            // Execute it with the provided parameters
            results result;
            conn.execute_statement(stmt, std::make_tuple(), result);

            // Validate the received contents
            validate_meta(result.meta(), table.metas);
            table.validate_rows(result.rows());
        }
    }
}

BOOST_FIXTURE_TEST_CASE(statement_write, database_types_fixture)
{
    start_transaction();

    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table.name)
        {
            // Prepare the statements
            auto insert_stmt = conn.prepare_statement(table.insert_sql());
            auto query_stmt = conn.prepare_statement(table.select_sql());

            // Remove all contents from the table
            results result;
            conn.query(table.delete_sql(), result);

            // Insert all the contents again
            boost::mysql::execution_state st;
            for (const auto& row : table.rws)
            {
                conn.start_statement_execution(insert_stmt, row.begin(), row.end(), st);
                BOOST_TEST_REQUIRE(st.complete());
            }

            // Query them again and verify the insertion was okay
            conn.execute_statement(query_stmt, std::make_tuple(), result);
            validate_meta(result.meta(), table.metas);
            table.validate_rows(result.rows());
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()  // test_database_types
