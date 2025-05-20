//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
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
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/column_type.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows_view.hpp>
#include <boost/mysql/sequence.hpp>
#include <boost/mysql/static_results.hpp>
#include <boost/mysql/tcp.hpp>

#include <boost/describe/class.hpp>
#include <boost/describe/members.hpp>
#include <boost/describe/modifiers.hpp>
#include <boost/describe/operators.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/none.hpp>
#include <boost/none_t.hpp>
#include <boost/optional/optional.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <ostream>
#include <stdint.h>
#include <type_traits>
#include <vector>

#include "test_common/create_basic.hpp"
#include "test_common/network_result.hpp"
#include "test_common/printing.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/metadata_validator.hpp"
#include "test_integration/server_features.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::none;
using boost::optional;

BOOST_AUTO_TEST_SUITE(test_database_types)

#ifdef BOOST_MYSQL_CXX14

// operator<< doesn't work for blob (vector<unsigned char>) or time (chrono::duration<...>)
template <class T, class = typename std::enable_if<boost::describe::has_describe_members<T>::value>::type>
std::ostream& operator<<(std::ostream& os, const T& v)
{
    os << '{';
    boost::mp11::mp_for_each<boost::describe::describe_members<T, boost::describe::mod_public>>([&](auto D) {
        os << "." << D.name << " = " << boost::mysql::detail::to_field(v.*D.pointer) << ", ";
    });
    return os << '}';
}
using boost::describe::operators::operator==;

#endif

// Helpers
using flagsvec = std::vector<meta_validator::flag_getter>;

const flagsvec flags_unsigned{&metadata::is_unsigned};
const flagsvec flags_zerofill{&metadata::is_unsigned, &metadata::is_zerofill};
const flagsvec no_flags{};

constexpr format_options opts{utf8mb4_charset, true};

struct database_types_fixture : any_connection_fixture
{
    database_types_fixture()
    {
        // Connect
        connect(connect_params_builder().multi_queries(true).build());

        // Sets the time_zone to a well known value, so we can deterministically read TIMESTAMPs
        // Sets also sql_mode to allow invalid dates
        results result;
        conn.async_execute(
                "SET session time_zone = '+02:00'; "
                "SET session sql_mode = 'ALLOW_INVALID_DATES'",
                result,
                as_netresult
        )
            .validate_no_error();
    }
};

struct table_base
{
    std::string name;
    std::vector<meta_validator> metas;
    std::vector<boost::mysql::row> rws;

    table_base(std::string name) : name(std::move(name))
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

    virtual ~table_base() {}
    virtual void select_static(any_connection& conn) = 0;

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

    void validate_rows(rows_view actual) const
    {
        // Sort the expected rows as the database retrieves it
        std::vector<row> expected{rws.begin(), rws.end()};
        std::sort(expected.begin(), expected.end(), [](row_view r1, row_view r2) {
            return r1.at(0).as_string() < r2.at(0).as_string();
        });

        // Compare
        BOOST_TEST(actual == expected, boost::test_tools::per_element());
    }

    std::string select_sql() const { return format_sql(opts, "SELECT * FROM {:i} ORDER BY id", name); }

    std::string insert_sql_stmt() const
    {
        auto format_fn = [](const meta_validator&, format_context_base& ctx) { ctx.append_raw("?"); };
        return format_sql(opts, "INSERT INTO {:i} VALUES ({})", name, sequence(metas, format_fn));
    }

    std::string insert_sql() const
    {
        auto format_fn = [](row_view r, format_context_base& ctx) { format_sql_to(ctx, "({})", r); };
        return format_sql(opts, "INSERT INTO {:i} VALUES {}", name, sequence(rws, format_fn));
    }

    std::string delete_sql() const { return format_sql(opts, "DELETE FROM {:i}", name); }
};

template <class T>
T convert_none(T v) noexcept
{
    return v;
}

std::nullptr_t convert_none(boost::none_t) { return nullptr; }

template <class StaticRow>
class table : public table_base
{
    std::vector<StaticRow> static_rows_;

public:
    using table_base::table_base;

    template <typename... Args>
    void add_row(const char* id, const Args&... args)
    {
        assert(sizeof...(Args) + 1 == metas.size());
        static_rows_.push_back(StaticRow{id, args...});
        rws.emplace_back(makerow(id, convert_none(args)...));
    }

#ifdef BOOST_MYSQL_CXX14
    void select_static(any_connection& conn) override
    {
        // Execute the query
        static_results<StaticRow> result;
        conn.async_execute(select_sql(), result, as_netresult).validate_no_error();

        // Validate metadata
        validate_meta(result.meta(), metas);

        // Validate the rows
        std::vector<StaticRow> expected{static_rows_};
        std::sort(expected.begin(), expected.end(), [](const StaticRow& r1, const StaticRow& r2) {
            return r1.id < r2.id;
        });
        BOOST_TEST(result.rows() == expected, boost::test_tools::per_element());
    }
#else
    void select_static(any_connection&) override {}
#endif
};

template <class StaticRow>
std::unique_ptr<table<StaticRow>> make_table(std::string name)
{
    return std::unique_ptr<table<StaticRow>>(new table<StaticRow>(std::move(name)));
}

using table_ptr = std::unique_ptr<table_base>;

// Int helpers
void int_table_columns(table_base& output, column_type type)
{
    output.add_meta("field_signed", type);
    output.add_meta("field_unsigned", type, flags_unsigned);
    output.add_meta("field_width", type);
    output.add_meta("field_zerofill", type, flags_zerofill);
}

template <class SignedInt, class UnsignedInt>
struct int_table_row
{
    std::string id;
    optional<SignedInt> field_signed;
    optional<UnsignedInt> field_unsigned;
    optional<SignedInt> field_width;
    optional<UnsignedInt> field_zerofill;
};

// TINYINT
using tinyint_row = int_table_row<int8_t, uint8_t>;
BOOST_DESCRIBE_STRUCT(tinyint_row, (), (id, field_signed, field_unsigned, field_width, field_zerofill))

static table_ptr types_tinyint()
{
    auto res = make_table<tinyint_row>("types_tinyint");
    int_table_columns(*res, column_type::tinyint);

    // clang-format off
    res->add_row("regular",  int8_t(20),    uint8_t(20),    int8_t(20),  uint8_t(20));
    res->add_row("negative", int8_t(-20),   none,           int8_t(-20), none);
    res->add_row("min",      int8_t(-0x80), uint8_t(0u),    none,        uint8_t(0u));
    res->add_row("max",      int8_t(0x7f),  uint8_t(0xffu), none,        none);
    // clang-format on
    return table_ptr(std::move(res));
}

// SMALLINT
using smallint_row = int_table_row<int16_t, uint16_t>;
BOOST_DESCRIBE_STRUCT(smallint_row, (), (id, field_signed, field_unsigned, field_width, field_zerofill))

static table_ptr types_smallint()
{
    auto res = make_table<smallint_row>("types_smallint");
    int_table_columns(*res, column_type::smallint);

    // clang-format off
    res->add_row("regular",  int16_t(20),      uint16_t(20u),     int16_t(20),  uint16_t(20u));
    res->add_row("negative", int16_t(-20),     none,              int16_t(-20), none);
    res->add_row("min",      int16_t(-0x8000), uint16_t(0u),      none,         uint16_t(0u));
    res->add_row("max",      int16_t(0x7fff),  uint16_t(0xffffu), none,         none);
    // clang-format on
    return table_ptr(std::move(res));  // This cast is required due to a bug in old clangs
}

// MEDIUMINT (row type shared with INT)
using int_row = int_table_row<int32_t, uint32_t>;
BOOST_DESCRIBE_STRUCT(int_row, (), (id, field_signed, field_unsigned, field_width, field_zerofill))

static table_ptr types_mediumint()
{
    auto res = make_table<int_row>("types_mediumint");
    int_table_columns(*res, column_type::mediumint);

    // clang-format off
    res->add_row("regular",  int32_t(20),        uint32_t(20u),       int32_t(20),  uint32_t(20u));
    res->add_row("negative", int32_t(-20),       none,                int32_t(-20), none);
    res->add_row("min",      int32_t(-0x800000), uint32_t(0u),        none,         uint32_t(0u));
    res->add_row("max",      int32_t(0x7fffff),  uint32_t(0xffffffu), none,         none);
    // clang-format on
    return table_ptr(std::move(res));
}

// INT
static table_ptr types_int()
{
    auto res = make_table<int_row>("types_int");
    int_table_columns(*res, column_type::int_);

    // clang-format off
    res->add_row("regular",  int32_t(20),            uint32_t(20u),         int32_t(20),  uint32_t(20u));
    res->add_row("negative", int32_t(-20),           none,                  int32_t(-20), none);
    res->add_row("min",      int32_t(-0x80000000LL), uint32_t(0u),          none,         0u);
    res->add_row("max",      int32_t(0x7fffffff),    uint32_t(0xffffffffu), none,         none);
    // clang-format on
    return table_ptr(std::move(res));
}

// BIGINT
using bigint_row = int_table_row<int64_t, uint64_t>;
BOOST_DESCRIBE_STRUCT(bigint_row, (), (id, field_signed, field_unsigned, field_width, field_zerofill))

static table_ptr types_bigint()
{
    auto res = make_table<bigint_row>("types_bigint");
    int_table_columns(*res, column_type::bigint);

    // clang-format off
    res->add_row("regular",   20,                     20u,                 20,   20u);
    res->add_row("negative", -20,                     none,                -20,  none);
    res->add_row("min",      -0x7fffffffffffffff - 1, 0u,                  none, 0u);
    res->add_row("max",       0x7fffffffffffffff,     0xffffffffffffffffu, none, none);
    // clang-format on
    return table_ptr(std::move(res));
}

// YEAR
struct year_row
{
    std::string id;
    optional<uint16_t> field_default;
};
BOOST_DESCRIBE_STRUCT(year_row, (), (id, field_default))

static table_ptr types_year()
{
    auto res = make_table<year_row>("types_year");
    res->add_meta("field_default", column_type::year, flags_zerofill);

    // clang-format off
    res->add_row("regular", uint16_t(2019u));
    res->add_row("min",     uint16_t(1901u));
    res->add_row("max",     uint16_t(2155u));
    res->add_row("zero",    uint16_t(0u));
    // clang-format on
    return table_ptr(std::move(res));
}

// BOOL
struct bool_row
{
    std::string id;
    optional<bool> field_default;
};
BOOST_DESCRIBE_STRUCT(bool_row, (), (id, field_default))

static table_ptr types_bool()
{
    auto res = make_table<bool_row>("types_bool");
    res->add_meta("field_default", column_type::tinyint);

    res->add_row("true", true);
    res->add_row("false", false);
    return table_ptr(std::move(res));
}

// BIT
struct bit_row
{
    std::string id;
    optional<uint64_t> field_1;
    optional<uint64_t> field_8;
    optional<uint64_t> field_14;
    optional<uint64_t> field_16;
    optional<uint64_t> field_24;
    optional<uint64_t> field_25;
    optional<uint64_t> field_32;
    optional<uint64_t> field_40;
    optional<uint64_t> field_48;
    optional<uint64_t> field_56;
    optional<uint64_t> field_64;
};
BOOST_DESCRIBE_STRUCT(
    bit_row,
    (),
    (id,
     field_1,
     field_8,
     field_14,
     field_16,
     field_24,
     field_25,
     field_32,
     field_40,
     field_48,
     field_56,
     field_64)
)

static table_ptr types_bit()
{
    auto res = make_table<bit_row>("types_bit");
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
        "field_64"
    };
    for (const char* col : columns)
        res->add_meta(col, column_type::bit, flags_unsigned);

    // clang-format off
    res->add_row("min",     0u, 0u,    0u,      0u,      0u,        0u,         0u,          0u,            0u,              0u,                0u);
    res->add_row("regular", 1u, 0x9eu, 0x1e2au, 0x1234u, 0x123456u, 0x154abe0u, 0x12345678u, 0x123456789au, 0x123456789abcu, 0x123456789abcdeu, 0x1234567812345678u);
    res->add_row("max",     1u, 0xffu, 0x3fffu, 0xffffu, 0xffffffu, 0x1ffffffu, 0xffffffffu, 0xffffffffffu, 0xffffffffffffu, 0xffffffffffffffu, 0xffffffffffffffffu);
    // clang-format on
    return table_ptr(std::move(res));
}

// FLOAT
struct float_row
{
    std::string id;
    optional<float> field_signed;
    optional<float> field_unsigned;
    optional<float> field_width;
    optional<float> field_zerofill;
};
BOOST_DESCRIBE_STRUCT(float_row, (), (id, field_signed, field_unsigned, field_width, field_zerofill))

static table_ptr types_float()
{
    auto res = make_table<float_row>("types_float");
    res->add_meta("field_signed", column_type::float_, no_flags, 31);
    res->add_meta("field_unsigned", column_type::float_, flags_unsigned, 31);
    res->add_meta("field_width", column_type::float_, no_flags, 10);
    res->add_meta("field_zerofill", column_type::float_, flags_zerofill, 31);

    // clang-format off
    res->add_row("zero",                              0.f,       0.f,   0.f,  0.f);
    res->add_row("int_positive",                      4.f,       none,  none, none);
    res->add_row("int_negative",                     -4.f,       none,  none, none);
    res->add_row("fractional_positive",               4.2f,      4.2f,  4.2f, 4.2f);
    res->add_row("fractional_negative",              -4.2f,      none, -4.2f, none);
    res->add_row("positive_exp_positive_int",         3e20f,     none,  none, none);
    res->add_row("positive_exp_negative_int",        -3e20f,     none,  none, none);
    res->add_row("positive_exp_positive_fractional",  3.14e20f,  none,  none, 3.14e20f);
    res->add_row("positive_exp_negative_fractional", -3.14e20f,  none,  none, none);
    res->add_row("negative_exp_positive_fractional",  3.14e-20f, none,  none, 3.14e-20f);
    // clang-format on
    return table_ptr(std::move(res));
}

// DOUBLE
struct double_row
{
    std::string id;
    optional<double> field_signed;
    optional<double> field_unsigned;
    optional<double> field_width;
    optional<double> field_zerofill;
};
BOOST_DESCRIBE_STRUCT(double_row, (), (id, field_signed, field_unsigned, field_width, field_zerofill))

static table_ptr types_double()
{
    auto res = make_table<double_row>("types_double");
    res->add_meta("field_signed", column_type::double_, no_flags, 31);
    res->add_meta("field_unsigned", column_type::double_, flags_unsigned, 31);
    res->add_meta("field_width", column_type::double_, no_flags, 10);
    res->add_meta("field_zerofill", column_type::double_, flags_zerofill, 31);

    // clang-format off
    res->add_row("zero",                              0.,        0.,   0.,   0.);
    res->add_row("int_positive",                      4.,        none, none, none);
    res->add_row("int_negative",                     -4.,        none, none, none);
    res->add_row("fractional_positive",               4.2,       4.2,   4.2, 4.2);
    res->add_row("fractional_negative",              -4.2,       none, -4.2, none);
    res->add_row("positive_exp_positive_int",         3e200,     none, none, none);
    res->add_row("positive_exp_negative_int",        -3e200,     none, none, none);
    res->add_row("positive_exp_positive_fractional",  3.14e200,  none, none, 3.14e200);
    res->add_row("positive_exp_negative_fractional", -3.14e200,  none, none, none);
    res->add_row("negative_exp_positive_fractional",  3.14e-200, none, none, 3.14e-200);
    // clang-format on
    return table_ptr(std::move(res));
}

// DATE
struct date_row
{
    std::string id;
    optional<date> field_date;
};
BOOST_DESCRIBE_STRUCT(date_row, (), (id, field_date))

static table_ptr types_date()
{
    auto res = make_table<date_row>("types_date");
    res->add_meta("field_date", column_type::date);

    // clang-format off
    res->add_row("regular",                           date(2010u, 3u, 28u));
    res->add_row("leap_regular",                      date(1788u, 2u, 29u));
    res->add_row("leap_400",                          date(2000u, 2u, 29u));
    res->add_row("min",                               date(0u, 1u, 01u));
    res->add_row("max",                               date(9999u, 12u, 31u));
    res->add_row("zero",                              date(0u, 0u, 0u));
    res->add_row("yzero_mzero_dregular",              date(0u, 0u, 20u));
    res->add_row("yzero_mregular_dzero",              date(0u, 11u, 0u));
    res->add_row("yzero_invalid_date",                date(0u, 11u, 31u));
    res->add_row("yregular_mzero_dzero",              date(2020u, 0u, 0u));
    res->add_row("yregular_mzero_dregular",           date(2020u, 0u, 20u));
    res->add_row("yregular_mregular_dzero",           date(2020u, 11u, 0u));
    res->add_row("yregular_invalid_date",             date(2020u, 11u, 31u));
    res->add_row("yregular_invalid_date_leapregular", date(1999u, 2u, 29u));
    res->add_row("yregular_invalid_date_leap100",     date(1900u, 2u, 29u));
    // clang-format on
    return table_ptr(std::move(res));
}

// DATETIME and TIMESTAMP (they share row definition)
struct datetime_row
{
    std::string id;
    optional<datetime> field_0;
    optional<datetime> field_1;
    optional<datetime> field_2;
    optional<datetime> field_3;
    optional<datetime> field_4;
    optional<datetime> field_5;
    optional<datetime> field_6;
};
BOOST_DESCRIBE_STRUCT(datetime_row, (), (id, field_0, field_1, field_2, field_3, field_4, field_5, field_6))

void datetime_timestamp_common_rows(table<datetime_row>& res)
{
    // clang-format off
    res.add_row("date",         datetime(2010, 5, 2, 0, 0, 0, 0),   datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0),           datetime(2010, 5, 2, 0, 0, 0, 0));
    res.add_row("date_leap4",   datetime(2004, 2, 29, 0, 0, 0, 0),  datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0),          datetime(2004, 2, 29, 0, 0, 0, 0));
    res.add_row("date_leap400", datetime(2000, 2, 29, 0, 0, 0, 0),  datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0),          datetime(2000, 2, 29, 0, 0, 0, 0));
    res.add_row("u",            none,                               datetime(2010, 5, 2, 0, 0, 0, 100000),      datetime(2010, 5, 2, 0, 0, 0, 120000),      datetime(2010, 5, 2, 0, 0, 0, 123000),      datetime(2010, 5, 2, 0, 0, 0, 123400),      datetime(2010, 5, 2, 0, 0, 0, 123450),      datetime(2010, 5, 2, 0, 0, 0, 123456));
    res.add_row("s",            datetime(2010, 5, 2, 0, 0, 50, 0),  datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0),          datetime(2010, 5, 2, 0, 0, 50, 0)), 
    res.add_row("m",            datetime(2010, 5, 2, 0, 1, 0, 0),   datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0),           datetime(2010, 5, 2, 0, 1, 0, 0));
    res.add_row("hs",           datetime(2010, 5, 2, 23, 0, 50, 0), datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0),         datetime(2010, 5, 2, 23, 0, 50, 0));
    res.add_row("ms",           datetime(2010, 5, 2, 0, 1, 50, 0),  datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0),          datetime(2010, 5, 2, 0, 1, 50, 0));
    res.add_row("hu",           none,                               datetime(2010, 5, 2, 23, 0, 0, 100000),     datetime(2010, 5, 2, 23, 0, 0, 120000),     datetime(2010, 5, 2, 23, 0, 0, 123000),     datetime(2010, 5, 2, 23, 0, 0, 123400),     datetime(2010, 5, 2, 23, 0, 0, 123450),     datetime(2010, 5, 2, 23, 0, 0, 123456));
    res.add_row("mu",           none,                               datetime(2010, 5, 2, 0, 1, 0, 100000),      datetime(2010, 5, 2, 0, 1, 0, 120000),      datetime(2010, 5, 2, 0, 1, 0, 123000),      datetime(2010, 5, 2, 0, 1, 0, 123400),      datetime(2010, 5, 2, 0, 1, 0, 123450),      datetime(2010, 5, 2, 0, 1, 0, 123456));
    res.add_row("hmu",          none,                               datetime(2010, 5, 2, 23, 1, 0, 100000),     datetime(2010, 5, 2, 23, 1, 0, 120000),     datetime(2010, 5, 2, 23, 1, 0, 123000),     datetime(2010, 5, 2, 23, 1, 0, 123400),     datetime(2010, 5, 2, 23, 1, 0, 123450),     datetime(2010, 5, 2, 23, 1, 0, 123456));
    res.add_row("su",           none,                               datetime(2010, 5, 2, 0, 0, 50, 100000),     datetime(2010, 5, 2, 0, 0, 50, 120000),     datetime(2010, 5, 2, 0, 0, 50, 123000),     datetime(2010, 5, 2, 0, 0, 50, 123400),     datetime(2010, 5, 2, 0, 0, 50, 123450),     datetime(2010, 5, 2, 0, 0, 50, 123456));
    res.add_row("hsu",          none,                               datetime(2010, 5, 2, 23, 0, 50, 100000),    datetime(2010, 5, 2, 23, 0, 50, 120000),    datetime(2010, 5, 2, 23, 0, 50, 123000),    datetime(2010, 5, 2, 23, 0, 50, 123400),    datetime(2010, 5, 2, 23, 0, 50, 123450),    datetime(2010, 5, 2, 23, 0, 50, 123456));
    res.add_row("msu",          none,                               datetime(2010, 5, 2, 0, 1, 50, 100000),     datetime(2010, 5, 2, 0, 1, 50, 120000),     datetime(2010, 5, 2, 0, 1, 50, 123000),     datetime(2010, 5, 2, 0, 1, 50, 123400),     datetime(2010, 5, 2, 0, 1, 50, 123450),     datetime(2010, 5, 2, 0, 1, 50, 123456));
    res.add_row("h",            datetime(2010, 5, 2, 23, 0, 0, 0),  datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0),          datetime(2010, 5, 2, 23, 0, 0, 0));
    res.add_row("hm",           datetime(2010, 5, 2, 23, 1, 0, 0),  datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0),          datetime(2010, 5, 2, 23, 1, 0, 0));
    res.add_row("hms",          datetime(2010, 5, 2, 23, 1, 50, 0), datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0),         datetime(2010, 5, 2, 23, 1, 50, 0));
    res.add_row("hmsu",         none,                               datetime(2010, 5, 2, 23, 1, 50, 100000),    datetime(2010, 5, 2, 23, 1, 50, 120000),    datetime(2010, 5, 2, 23, 1, 50, 123000),    datetime(2010, 5, 2, 23, 1, 50, 123400),    datetime(2010, 5, 2, 23, 1, 50, 123450),    datetime(2010, 5, 2, 23, 1, 50, 123456));
    // clang-format on
}

static table_ptr types_datetime()
{
    auto res = make_table<datetime_row>("types_datetime");
    res->add_meta("field_0", column_type::datetime, no_flags, 0, flags_unsigned);
    res->add_meta("field_1", column_type::datetime, no_flags, 1, flags_unsigned);
    res->add_meta("field_2", column_type::datetime, no_flags, 2, flags_unsigned);
    res->add_meta("field_3", column_type::datetime, no_flags, 3, flags_unsigned);
    res->add_meta("field_4", column_type::datetime, no_flags, 4, flags_unsigned);
    res->add_meta("field_5", column_type::datetime, no_flags, 5, flags_unsigned);
    res->add_meta("field_6", column_type::datetime, no_flags, 6, flags_unsigned);

    datetime_timestamp_common_rows(*res);

    // clang-format off
    res->add_row("min",                                     datetime(0, 1, 1, 0, 0, 0, 0),         datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0),              datetime(0, 1, 1, 0, 0, 0, 0));
    res->add_row("max",                                     datetime(9999, 12, 31, 23, 59, 59, 0), datetime(9999, 12, 31, 23, 59, 59, 900000), datetime(9999, 12, 31, 23, 59, 59, 990000), datetime(9999, 12, 31, 23, 59, 59, 999000), datetime(9999, 12, 31, 23, 59, 59, 999900), datetime(9999, 12, 31, 23, 59, 59, 999990), datetime(9999, 12, 31, 23, 59, 59, 999999));
    res->add_row("date_zero",                               datetime(   0,  0,  0,  0,  0,  0, 0), datetime(   0,  0,  0,  0,  0,  0, 0),      datetime(   0,  0,  0,  0,  0,  0,  0),     datetime(   0,  0,  0,  0,  0,  0,   0),    datetime(   0,  0,  0,  0,  0,  0,    0),   datetime(   0,  0,  0,  0,  0,  0,     0),  datetime(   0,  0,  0,  0,  0,  0,      0));
    res->add_row("date_yzero_mzero_dregular",               datetime(   0,  0, 10,  0,  0,  0, 0), datetime(   0,  0, 10,  0,  0,  0, 0),      datetime(   0,  0, 10,  0,  0,  0,  0),     datetime(   0,  0, 10,  0,  0,  0,   0),    datetime(   0,  0, 10,  0,  0,  0,    0),   datetime(   0,  0, 10,  0,  0,  0,     0),  datetime(   0,  0, 10,  0,  0,  0,      0));
    res->add_row("date_yzero_mregular_dzero",               datetime(   0, 10,  0,  0,  0,  0, 0), datetime(   0, 10,  0,  0,  0,  0, 0),      datetime(   0, 10,  0,  0,  0,  0,  0),     datetime(   0, 10,  0,  0,  0,  0,   0),    datetime(   0, 10,  0,  0,  0,  0,    0),   datetime(   0, 10,  0,  0,  0,  0,     0),  datetime(   0, 10,  0,  0,  0,  0,      0));
    res->add_row("date_yzero_invalid_date",                 datetime(   0, 11, 31,  0,  0,  0, 0), datetime(   0, 11, 31,  0,  0,  0, 0),      datetime(   0, 11, 31,  0,  0,  0,  0),     datetime(   0, 11, 31,  0,  0,  0,   0),    datetime(   0, 11, 31,  0,  0,  0,    0),   datetime(   0, 11, 31,  0,  0,  0,     0),  datetime(   0, 11, 31,  0,  0,  0,      0));
    res->add_row("date_yregular_mzero_dzero",               datetime(2020,  0,  0,  0,  0,  0, 0), datetime(2020,  0,  0,  0,  0,  0, 0),      datetime(2020,  0,  0,  0,  0,  0,  0),     datetime(2020,  0,  0,  0,  0,  0,   0),    datetime(2020,  0,  0,  0,  0,  0,    0),   datetime(2020,  0,  0,  0,  0,  0,     0),  datetime(2020,  0,  0,  0,  0,  0,      0));
    res->add_row("date_yregular_mzero_dregular",            datetime(2020,  0, 10,  0,  0,  0, 0), datetime(2020,  0, 10,  0,  0,  0, 0),      datetime(2020,  0, 10,  0,  0,  0,  0),     datetime(2020,  0, 10,  0,  0,  0,   0),    datetime(2020,  0, 10,  0,  0,  0,    0),   datetime(2020,  0, 10,  0,  0,  0,     0),  datetime(2020,  0, 10,  0,  0,  0,      0));
    res->add_row("date_yregular_mregular_dzero",            datetime(2020, 10,  0,  0,  0,  0, 0), datetime(2020, 10,  0,  0,  0,  0, 0),      datetime(2020, 10,  0,  0,  0,  0,  0),     datetime(2020, 10,  0,  0,  0,  0,   0),    datetime(2020, 10,  0,  0,  0,  0,    0),   datetime(2020, 10,  0,  0,  0,  0,     0),  datetime(2020, 10,  0,  0,  0,  0,      0));
    res->add_row("date_yregular_invalid_date",              datetime(2020, 11, 31,  0,  0,  0, 0), datetime(2020, 11, 31,  0,  0,  0, 0),      datetime(2020, 11, 31,  0,  0,  0,  0),     datetime(2020, 11, 31,  0,  0,  0,   0),    datetime(2020, 11, 31,  0,  0,  0,    0),   datetime(2020, 11, 31,  0,  0,  0,     0),  datetime(2020, 11, 31,  0,  0,  0,      0));
    res->add_row("date_yregular_invalid_date_leapregular",  datetime(1999,  2, 29,  0,  0,  0, 0), datetime(1999,  2, 29,  0,  0,  0, 0),      datetime(1999,  2, 29,  0,  0,  0,  0),     datetime(1999,  2, 29,  0,  0,  0,   0),    datetime(1999,  2, 29,  0,  0,  0,    0),   datetime(1999,  2, 29,  0,  0,  0,     0),  datetime(1999,  2, 29,  0,  0,  0,      0));
    res->add_row("date_yregular_invalid_date_leap100",      datetime(1900,  2, 29,  0,  0,  0, 0), datetime(1900,  2, 29,  0,  0,  0, 0),      datetime(1900,  2, 29,  0,  0,  0,  0),     datetime(1900,  2, 29,  0,  0,  0,   0),    datetime(1900,  2, 29,  0,  0,  0,    0),   datetime(1900,  2, 29,  0,  0,  0,     0),  datetime(1900,  2, 29,  0,  0,  0,      0)),
    res->add_row("hms_zero",                                datetime(   0,  0,  0, 10, 20, 30, 0), datetime(   0,  0,  0, 10, 20, 30, 0),      datetime(   0,  0,  0, 10, 20, 30,  0),     datetime(   0,  0,  0, 10, 20, 30,   0),    datetime(   0,  0,  0, 10, 20, 30,    0),   datetime(   0,  0,  0, 10, 20, 30,     0),  datetime(   0,  0,  0, 10, 20, 30,      0));
    res->add_row("hms_yzero_mzero_dregular",                datetime(   0,  0, 10, 10, 20, 30, 0), datetime(   0,  0, 10, 10, 20, 30, 0),      datetime(   0,  0, 10, 10, 20, 30,  0),     datetime(   0,  0, 10, 10, 20, 30,   0),    datetime(   0,  0, 10, 10, 20, 30,    0),   datetime(   0,  0, 10, 10, 20, 30,     0),  datetime(   0,  0, 10, 10, 20, 30,      0));
    res->add_row("hms_yzero_mregular_dzero",                datetime(   0, 10,  0, 10, 20, 30, 0), datetime(   0, 10,  0, 10, 20, 30, 0),      datetime(   0, 10,  0, 10, 20, 30,  0),     datetime(   0, 10,  0, 10, 20, 30,   0),    datetime(   0, 10,  0, 10, 20, 30,    0),   datetime(   0, 10,  0, 10, 20, 30,     0),  datetime(   0, 10,  0, 10, 20, 30,      0));
    res->add_row("hms_yzero_invalid_date",                  datetime(   0, 11, 31, 10, 20, 30, 0), datetime(   0, 11, 31, 10, 20, 30, 0),      datetime(   0, 11, 31, 10, 20, 30,  0),     datetime(   0, 11, 31, 10, 20, 30,   0),    datetime(   0, 11, 31, 10, 20, 30,    0),   datetime(   0, 11, 31, 10, 20, 30,     0),  datetime(   0, 11, 31, 10, 20, 30,      0));
    res->add_row("hms_yregular_mzero_dzero",                datetime(2020,  0,  0, 10, 20, 30, 0), datetime(2020,  0,  0, 10, 20, 30, 0),      datetime(2020,  0,  0, 10, 20, 30,  0),     datetime(2020,  0,  0, 10, 20, 30,   0),    datetime(2020,  0,  0, 10, 20, 30,    0),   datetime(2020,  0,  0, 10, 20, 30,     0),  datetime(2020,  0,  0, 10, 20, 30,      0));
    res->add_row("hms_yregular_mzero_dregular",             datetime(2020,  0, 10, 10, 20, 30, 0), datetime(2020,  0, 10, 10, 20, 30, 0),      datetime(2020,  0, 10, 10, 20, 30,  0),     datetime(2020,  0, 10, 10, 20, 30,   0),    datetime(2020,  0, 10, 10, 20, 30,    0),   datetime(2020,  0, 10, 10, 20, 30,     0),  datetime(2020,  0, 10, 10, 20, 30,      0));
    res->add_row("hms_yregular_mregular_dzero",             datetime(2020, 10,  0, 10, 20, 30, 0), datetime(2020, 10,  0, 10, 20, 30, 0),      datetime(2020, 10,  0, 10, 20, 30,  0),     datetime(2020, 10,  0, 10, 20, 30,   0),    datetime(2020, 10,  0, 10, 20, 30,    0),   datetime(2020, 10,  0, 10, 20, 30,     0),  datetime(2020, 10,  0, 10, 20, 30,      0));
    res->add_row("hms_yregular_invalid_date",               datetime(2020, 11, 31, 10, 20, 30, 0), datetime(2020, 11, 31, 10, 20, 30, 0),      datetime(2020, 11, 31, 10, 20, 30,  0),     datetime(2020, 11, 31, 10, 20, 30,   0),    datetime(2020, 11, 31, 10, 20, 30,    0),   datetime(2020, 11, 31, 10, 20, 30,     0),  datetime(2020, 11, 31, 10, 20, 30,      0));
    res->add_row("hms_yregular_invalid_date_leapregular",   datetime(1999,  2, 29, 10, 20, 30, 0), datetime(1999,  2, 29, 10, 20, 30, 0),      datetime(1999,  2, 29, 10, 20, 30,  0),     datetime(1999,  2, 29, 10, 20, 30,   0),    datetime(1999,  2, 29, 10, 20, 30,    0),   datetime(1999,  2, 29, 10, 20, 30,     0),  datetime(1999,  2, 29, 10, 20, 30,      0));
    res->add_row("hms_yregular_invalid_date_leap100",       datetime(1900,  2, 29, 10, 20, 30, 0), datetime(1900,  2, 29, 10, 20, 30, 0),      datetime(1900,  2, 29, 10, 20, 30,  0),     datetime(1900,  2, 29, 10, 20, 30,   0),    datetime(1900,  2, 29, 10, 20, 30,    0),   datetime(1900,  2, 29, 10, 20, 30,     0),  datetime(1900,  2, 29, 10, 20, 30,      0)),
    res->add_row("hmsu_zero",                               datetime(   0,  0,  0, 10, 20, 30, 0), datetime(   0,  0,  0, 10, 20, 30, 900000), datetime(   0,  0,  0, 10, 20, 30, 990000), datetime(   0,  0,  0, 10, 20, 30, 999000), datetime(   0,  0,  0, 10, 20, 30, 999900), datetime(   0,  0,  0, 10, 20, 30, 999990), datetime(   0,  0,  0, 10, 20, 30, 999999));
    res->add_row("hmsu_yzero_mzero_dregular",               datetime(   0,  0, 10, 10, 20, 30, 0), datetime(   0,  0, 10, 10, 20, 30, 900000), datetime(   0,  0, 10, 10, 20, 30, 990000), datetime(   0,  0, 10, 10, 20, 30, 999000), datetime(   0,  0, 10, 10, 20, 30, 999900), datetime(   0,  0, 10, 10, 20, 30, 999990), datetime(   0,  0, 10, 10, 20, 30, 999999));
    res->add_row("hmsu_yzero_mregular_dzero",               datetime(   0, 10,  0, 10, 20, 30, 0), datetime(   0, 10,  0, 10, 20, 30, 900000), datetime(   0, 10,  0, 10, 20, 30, 990000), datetime(   0, 10,  0, 10, 20, 30, 999000), datetime(   0, 10,  0, 10, 20, 30, 999900), datetime(   0, 10,  0, 10, 20, 30, 999990), datetime(   0, 10,  0, 10, 20, 30, 999999));
    res->add_row("hmsu_yzero_invalid_date",                 datetime(   0, 11, 31, 10, 20, 30, 0), datetime(   0, 11, 31, 10, 20, 30, 900000), datetime(   0, 11, 31, 10, 20, 30, 990000), datetime(   0, 11, 31, 10, 20, 30, 999000), datetime(   0, 11, 31, 10, 20, 30, 999900), datetime(   0, 11, 31, 10, 20, 30, 999990), datetime(   0, 11, 31, 10, 20, 30, 999999));
    res->add_row("hmsu_yregular_mzero_dzero",               datetime(2020,  0,  0, 10, 20, 30, 0), datetime(2020,  0,  0, 10, 20, 30, 900000), datetime(2020,  0,  0, 10, 20, 30, 990000), datetime(2020,  0,  0, 10, 20, 30, 999000), datetime(2020,  0,  0, 10, 20, 30, 999900), datetime(2020,  0,  0, 10, 20, 30, 999990), datetime(2020,  0,  0, 10, 20, 30, 999999));
    res->add_row("hmsu_yregular_mzero_dregular",            datetime(2020,  0, 10, 10, 20, 30, 0), datetime(2020,  0, 10, 10, 20, 30, 900000), datetime(2020,  0, 10, 10, 20, 30, 990000), datetime(2020,  0, 10, 10, 20, 30, 999000), datetime(2020,  0, 10, 10, 20, 30, 999900), datetime(2020,  0, 10, 10, 20, 30, 999990), datetime(2020,  0, 10, 10, 20, 30, 999999));
    res->add_row("hmsu_yregular_mregular_dzero",            datetime(2020, 10,  0, 10, 20, 30, 0), datetime(2020, 10,  0, 10, 20, 30, 900000), datetime(2020, 10,  0, 10, 20, 30, 990000), datetime(2020, 10,  0, 10, 20, 30, 999000), datetime(2020, 10,  0, 10, 20, 30, 999900), datetime(2020, 10,  0, 10, 20, 30, 999990), datetime(2020, 10,  0, 10, 20, 30, 999999));
    res->add_row("hmsu_yregular_invalid_date",              datetime(2020, 11, 31, 10, 20, 30, 0), datetime(2020, 11, 31, 10, 20, 30, 900000), datetime(2020, 11, 31, 10, 20, 30, 990000), datetime(2020, 11, 31, 10, 20, 30, 999000), datetime(2020, 11, 31, 10, 20, 30, 999900), datetime(2020, 11, 31, 10, 20, 30, 999990), datetime(2020, 11, 31, 10, 20, 30, 999999));
    res->add_row("hmsu_yregular_invalid_date_leapregular",  datetime(1999,  2, 29, 10, 20, 30, 0), datetime(1999,  2, 29, 10, 20, 30, 900000), datetime(1999,  2, 29, 10, 20, 30, 990000), datetime(1999,  2, 29, 10, 20, 30, 999000), datetime(1999,  2, 29, 10, 20, 30, 999900), datetime(1999,  2, 29, 10, 20, 30, 999990), datetime(1999,  2, 29, 10, 20, 30, 999999));
    res->add_row("hmsu_yregular_invalid_date_leap100",      datetime(1900,  2, 29, 10, 20, 30, 0), datetime(1900,  2, 29, 10, 20, 30, 900000), datetime(1900,  2, 29, 10, 20, 30, 990000), datetime(1900,  2, 29, 10, 20, 30, 999000), datetime(1900,  2, 29, 10, 20, 30, 999900), datetime(1900,  2, 29, 10, 20, 30, 999990), datetime(1900,  2, 29, 10, 20, 30, 999999));
    // clang-format on
    return table_ptr(std::move(res));
}

static table_ptr types_timestamp()
{
    auto res = make_table<datetime_row>("types_timestamp");
    res->add_meta("field_0", column_type::timestamp, no_flags, 0, flags_unsigned);
    res->add_meta("field_1", column_type::timestamp, no_flags, 1, flags_unsigned);
    res->add_meta("field_2", column_type::timestamp, no_flags, 2, flags_unsigned);
    res->add_meta("field_3", column_type::timestamp, no_flags, 3, flags_unsigned);
    res->add_meta("field_4", column_type::timestamp, no_flags, 4, flags_unsigned);
    res->add_meta("field_5", column_type::timestamp, no_flags, 5, flags_unsigned);
    res->add_meta("field_6", column_type::timestamp, no_flags, 6, flags_unsigned);

    datetime_timestamp_common_rows(*res);

    // clang-format off
    res->add_row("zero", datetime(),                            datetime(),                                 datetime(),                                 datetime(),                                 datetime(),                                 datetime(),                                 datetime());
    res->add_row("min",  datetime(1970,  1,  1,  2,  0,  1, 0), datetime(1970,  1,  1,  2,  0,  1, 0),      datetime(1970,  1,  1,  2,  0,  1,  0),     datetime(1970,  1,  1,  2,  0,  1,   0),    datetime(1970,  1,  1,  2,  0,  1,    0),   datetime(1970,  1,  1,  2,  0,  1,     0),  datetime(1970,  1,  1,  2,  0,  1,      0));
    res->add_row("max",  datetime(2038,  1, 19,  5, 14,  7, 0), datetime(2038,  1, 19,  5, 14,  7, 900000), datetime(2038,  1, 19,  5, 14,  7, 990000), datetime(2038,  1, 19,  5, 14,  7, 999000), datetime(2038,  1, 19,  5, 14,  7, 999900), datetime(2038,  1, 19,  5, 14,  7, 999990), datetime(2038,  1, 19,  5, 14,  7, 999999));
    // clang-format on
    return table_ptr(std::move(res));
}

// TIME
struct time_row
{
    std::string id;
    optional<boost::mysql::time> field_0;
    optional<boost::mysql::time> field_1;
    optional<boost::mysql::time> field_2;
    optional<boost::mysql::time> field_3;
    optional<boost::mysql::time> field_4;
    optional<boost::mysql::time> field_5;
    optional<boost::mysql::time> field_6;
};
BOOST_DESCRIBE_STRUCT(time_row, (), (id, field_0, field_1, field_2, field_3, field_4, field_5, field_6))

static table_ptr types_time()
{
    auto res = make_table<time_row>("types_time");
    res->add_meta("field_0", column_type::time, no_flags, 0, flags_unsigned);
    res->add_meta("field_1", column_type::time, no_flags, 1, flags_unsigned);
    res->add_meta("field_2", column_type::time, no_flags, 2, flags_unsigned);
    res->add_meta("field_3", column_type::time, no_flags, 3, flags_unsigned);
    res->add_meta("field_4", column_type::time, no_flags, 4, flags_unsigned);
    res->add_meta("field_5", column_type::time, no_flags, 5, flags_unsigned);
    res->add_meta("field_6", column_type::time, no_flags, 6, flags_unsigned);

    // clang-format off
    res->add_row("zero",            maket( 0,  0,  0),  maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0),          maket( 0,  0,  0));
    res->add_row("d",               maket(48,  0,  0),  maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0),          maket(48,  0,  0));
    res->add_row("negative_d",     -maket(48,  0,  0), -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0),         -maket(48,  0,  0));
    res->add_row("h",               maket(23,  0,  0),  maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0),          maket(23,  0,  0));
    res->add_row("negative_h",     -maket(23,  0,  0), -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0),         -maket(23,  0,  0));
    res->add_row("dh",              maket(71,  0,  0),  maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0),          maket(71,  0,  0));
    res->add_row("negative_dh",    -maket(71,  0,  0), -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0),         -maket(71,  0,  0));
    res->add_row("m",               maket( 0,  1,  0),  maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0),          maket( 0,  1,  0));
    res->add_row("negative_m",     -maket( 0,  1,  0), -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0),         -maket( 0,  1,  0));
    res->add_row("dm",              maket(48,  1,  0),  maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0),          maket(48,  1,  0));
    res->add_row("negative_dm",    -maket(48,  1,  0), -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0),         -maket(48,  1,  0));
    res->add_row("hm",              maket(23,  1,  0),  maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0),          maket(23,  1,  0));
    res->add_row("negative_hm",    -maket(23,  1,  0), -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0),         -maket(23,  1,  0));
    res->add_row("dhm",             maket(71,  1,  0),  maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0),          maket(71,  1,  0));
    res->add_row("negative_dhm",   -maket(71,  1,  0), -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0),         -maket(71,  1,  0));
    res->add_row("s",               maket( 0,  0, 50),  maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50),          maket( 0,  0, 50));
    res->add_row("negative_s",     -maket( 0,  0, 50), -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50),         -maket( 0,  0, 50));
    res->add_row("ds",              maket(48,  0, 50),  maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50),          maket(48,  0, 50));
    res->add_row("negative_ds",    -maket(48,  0, 50), -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50),         -maket(48,  0, 50));
    res->add_row("hs",              maket(23,  0, 50),  maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50),          maket(23,  0, 50));
    res->add_row("negative_hs",    -maket(23,  0, 50), -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50),         -maket(23,  0, 50));
    res->add_row("dhs",             maket(71,  0, 50),  maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50),          maket(71,  0, 50));
    res->add_row("negative_dhs",   -maket(71,  0, 50), -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50),         -maket(71,  0, 50));
    res->add_row("ms",              maket( 0,  1, 50),  maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50),          maket( 0,  1, 50));
    res->add_row("negative_ms",    -maket( 0,  1, 50), -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50),         -maket( 0,  1, 50));
    res->add_row("dms",             maket(48,  1, 50),  maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50),          maket(48,  1, 50));
    res->add_row("negative_dms",   -maket(48,  1, 50), -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50),         -maket(48,  1, 50));
    res->add_row("hms",             maket(23,  1, 50),  maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50),          maket(23,  1, 50));
    res->add_row("negative_hms",   -maket(23,  1, 50), -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50),         -maket(23,  1, 50));
    res->add_row("dhms",            maket(71,  1, 50),  maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50),          maket(71,  1, 50));
    res->add_row("negative_dhms",  -maket(71,  1, 50), -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50),         -maket(71,  1, 50));
    res->add_row("u",              none,                maket( 0,  0,  0, 100000),  maket( 0,  0,  0, 120000),  maket( 0,  0,  0, 123000),  maket( 0,  0,  0, 123400),  maket( 0,  0,  0, 123450),  maket( 0,  0,  0,  123456));
    res->add_row("negative_u",     none,               -maket( 0,  0,  0, 100000), -maket( 0,  0,  0, 120000), -maket( 0,  0,  0, 123000), -maket( 0,  0,  0, 123400), -maket( 0,  0,  0, 123450), -maket( 0,  0,  0,  123456));
    res->add_row("du",             none,                maket(48,  0,  0, 100000),  maket(48,  0,  0, 120000),  maket(48,  0,  0, 123000),  maket(48,  0,  0, 123400),  maket(48,  0,  0, 123450),  maket(48,  0,  0,  123456));
    res->add_row("negative_du",    none,               -maket(48,  0,  0, 100000), -maket(48,  0,  0, 120000), -maket(48,  0,  0, 123000), -maket(48,  0,  0, 123400), -maket(48,  0,  0, 123450), -maket(48,  0,  0,  123456));
    res->add_row("hu",             none,                maket(23,  0,  0, 100000),  maket(23,  0,  0, 120000),  maket(23,  0,  0, 123000),  maket(23,  0,  0, 123400),  maket(23,  0,  0, 123450),  maket(23,  0,  0,  123456));
    res->add_row("negative_hu",    none,               -maket(23,  0,  0, 100000), -maket(23,  0,  0, 120000), -maket(23,  0,  0, 123000), -maket(23,  0,  0, 123400), -maket(23,  0,  0, 123450), -maket(23,  0,  0,  123456));
    res->add_row("dhu",            none,                maket(71,  0,  0, 100000),  maket(71,  0,  0, 120000),  maket(71,  0,  0, 123000),  maket(71,  0,  0, 123400),  maket(71,  0,  0, 123450),  maket(71,  0,  0,  123456));
    res->add_row("negative_dhu",   none,               -maket(71,  0,  0, 100000), -maket(71,  0,  0, 120000), -maket(71,  0,  0, 123000), -maket(71,  0,  0, 123400), -maket(71,  0,  0, 123450), -maket(71,  0,  0,  123456));
    res->add_row("mu",             none,                maket( 0,  1,  0, 100000),  maket( 0,  1,  0, 120000),  maket( 0,  1,  0, 123000),  maket( 0,  1,  0, 123400),  maket( 0,  1,  0, 123450),  maket( 0,  1,  0,  123456));
    res->add_row("negative_mu",    none,               -maket( 0,  1,  0, 100000), -maket( 0,  1,  0, 120000), -maket( 0,  1,  0, 123000), -maket( 0,  1,  0, 123400), -maket( 0,  1,  0, 123450), -maket( 0,  1,  0,  123456));
    res->add_row("dmu",            none,                maket(48,  1,  0, 100000),  maket(48,  1,  0, 120000),  maket(48,  1,  0, 123000),  maket(48,  1,  0, 123400),  maket(48,  1,  0, 123450),  maket(48,  1,  0,  123456));
    res->add_row("negative_dmu",   none,               -maket(48,  1,  0, 100000), -maket(48,  1,  0, 120000), -maket(48,  1,  0, 123000), -maket(48,  1,  0, 123400), -maket(48,  1,  0, 123450), -maket(48,  1,  0,  123456));
    res->add_row("hmu",            none,                maket(23,  1,  0, 100000),  maket(23,  1,  0, 120000),  maket(23,  1,  0, 123000),  maket(23,  1,  0, 123400),  maket(23,  1,  0, 123450),  maket(23,  1,  0,  123456));
    res->add_row("negative_hmu",   none,               -maket(23,  1,  0, 100000), -maket(23,  1,  0, 120000), -maket(23,  1,  0, 123000), -maket(23,  1,  0, 123400), -maket(23,  1,  0, 123450), -maket(23,  1,  0,  123456));
    res->add_row("dhmu",           none,                maket(71,  1,  0, 100000),  maket(71,  1,  0, 120000),  maket(71,  1,  0, 123000),  maket(71,  1,  0, 123400),  maket(71,  1,  0, 123450),  maket(71,  1,  0,  123456));
    res->add_row("negative_dhmu",  none,               -maket(71,  1,  0, 100000), -maket(71,  1,  0, 120000), -maket(71,  1,  0, 123000), -maket(71,  1,  0, 123400), -maket(71,  1,  0, 123450), -maket(71,  1,  0,  123456));
    res->add_row("su",             none,                maket( 0,  0, 50, 100000),  maket( 0,  0, 50, 120000),  maket( 0,  0, 50, 123000),  maket( 0,  0, 50, 123400),  maket( 0,  0, 50, 123450),  maket( 0,  0, 50,  123456));
    res->add_row("negative_su",    none,               -maket( 0,  0, 50, 100000), -maket( 0,  0, 50, 120000), -maket( 0,  0, 50, 123000), -maket( 0,  0, 50, 123400), -maket( 0,  0, 50, 123450), -maket( 0,  0, 50,  123456));
    res->add_row("dsu",            none,                maket(48,  0, 50, 100000),  maket(48,  0, 50, 120000),  maket(48,  0, 50, 123000),  maket(48,  0, 50, 123400),  maket(48,  0, 50, 123450),  maket(48,  0, 50,  123456));
    res->add_row("negative_dsu",   none,               -maket(48,  0, 50, 100000), -maket(48,  0, 50, 120000), -maket(48,  0, 50, 123000), -maket(48,  0, 50, 123400), -maket(48,  0, 50, 123450), -maket(48,  0, 50,  123456));
    res->add_row("hsu",            none,                maket(23,  0, 50, 100000),  maket(23,  0, 50, 120000),  maket(23,  0, 50, 123000),  maket(23,  0, 50, 123400),  maket(23,  0, 50, 123450),  maket(23,  0, 50,  123456));
    res->add_row("negative_hsu",   none,               -maket(23,  0, 50, 100000), -maket(23,  0, 50, 120000), -maket(23,  0, 50, 123000), -maket(23,  0, 50, 123400), -maket(23,  0, 50, 123450), -maket(23,  0, 50,  123456));
    res->add_row("dhsu",           none,                maket(71,  0, 50, 100000),  maket(71,  0, 50, 120000),  maket(71,  0, 50, 123000),  maket(71,  0, 50, 123400),  maket(71,  0, 50, 123450),  maket(71,  0, 50,  123456));
    res->add_row("negative_dhsu",  none,               -maket(71,  0, 50, 100000), -maket(71,  0, 50, 120000), -maket(71,  0, 50, 123000), -maket(71,  0, 50, 123400), -maket(71,  0, 50, 123450), -maket(71,  0, 50,  123456));
    res->add_row("msu",            none,                maket( 0,  1, 50, 100000),  maket( 0,  1, 50, 120000),  maket( 0,  1, 50, 123000),  maket( 0,  1, 50, 123400),  maket( 0,  1, 50, 123450),  maket( 0,  1, 50,  123456));
    res->add_row("negative_msu",   none,               -maket( 0,  1, 50, 100000), -maket( 0,  1, 50, 120000), -maket( 0,  1, 50, 123000), -maket( 0,  1, 50, 123400), -maket( 0,  1, 50, 123450), -maket( 0,  1, 50,  123456));
    res->add_row("dmsu",           none,                maket(48,  1, 50, 100000),  maket(48,  1, 50, 120000),  maket(48,  1, 50, 123000),  maket(48,  1, 50, 123400),  maket(48,  1, 50, 123450),  maket(48,  1, 50,  123456));
    res->add_row("negative_dmsu",  none,               -maket(48,  1, 50, 100000), -maket(48,  1, 50, 120000), -maket(48,  1, 50, 123000), -maket(48,  1, 50, 123400), -maket(48,  1, 50, 123450), -maket(48,  1, 50,  123456));
    res->add_row("hmsu",           none,                maket(23,  1, 50, 100000),  maket(23,  1, 50, 120000),  maket(23,  1, 50, 123000),  maket(23,  1, 50, 123400),  maket(23,  1, 50, 123450),  maket(23,  1, 50,  123456));
    res->add_row("negative_hmsu",  none,               -maket(23,  1, 50, 100000), -maket(23,  1, 50, 120000), -maket(23,  1, 50, 123000), -maket(23,  1, 50, 123400), -maket(23,  1, 50, 123450), -maket(23,  1, 50,  123456));
    res->add_row("dhmsu",          none,                maket(71,  1, 50, 100000),  maket(71,  1, 50, 120000),  maket(71,  1, 50, 123000),  maket(71,  1, 50, 123400),  maket(71,  1, 50, 123450),  maket(71,  1, 50,  123456));
    res->add_row("negative_dhmsu", none,               -maket(71,  1, 50, 100000), -maket(71,  1, 50, 120000), -maket(71,  1, 50, 123000), -maket(71,  1, 50, 123400), -maket(71,  1, 50, 123450), -maket(71,  1, 50,  123456));
    res->add_row("min",            -maket(838, 59, 59),-maket(838,59, 58, 900000), -maket(838, 59, 58, 990000),-maket(838, 59, 58, 999000),-maket(838, 59, 58, 999900),-maket(838,59, 58, 999990), -maket(838,59, 58,  999999));
    res->add_row("max",             maket(838, 59, 59), maket(838,59, 58, 900000),  maket(838, 59, 58, 990000), maket(838, 59, 58, 999000), maket(838, 59, 58, 999900), maket(838,59, 58, 999990),  maket(838,59, 58,  999999));
    // clang-format on
    return table_ptr(std::move(res));
}

// string types
struct string_row
{
    std::string id;
    optional<std::string> field_char;
    optional<std::string> field_varchar;
    optional<std::string> field_tinytext;
    optional<std::string> field_text;
    optional<std::string> field_mediumtext;
    optional<std::string> field_longtext;
    optional<std::string> field_text_bincol;
    optional<std::string> field_enum;
    optional<std::string> field_set;
};
BOOST_DESCRIBE_STRUCT(
    string_row,
    (),
    (id,
     field_char,
     field_varchar,
     field_tinytext,
     field_text,
     field_mediumtext,
     field_longtext,
     field_text_bincol,
     field_enum,
     field_set)
)

static table_ptr types_string()
{
    auto res = make_table<string_row>("types_string");
    res->add_meta("field_char", column_type::char_);
    res->add_meta("field_varchar", column_type::varchar);
    res->add_meta("field_tinytext", column_type::text);
    res->add_meta("field_text", column_type::text);
    res->add_meta("field_mediumtext", column_type::text);
    res->add_meta("field_longtext", column_type::text);
    res->add_meta("field_text_bincol", column_type::text);
    res->add_meta("field_enum", column_type::enum_);
    res->add_meta("field_set", column_type::set);

    using std::string;

    // clang-format off
    res->add_row("regular", string("test_char"), string("test_varchar"), string("test_tinytext"), string("test_text"), string("test_mediumtext"), string("test_longtext"), string("test_bincol"), string("red"),    string("red,green"));
    res->add_row("utf8",    string("\xc3\xb1"),  string("\xc3\x91"),     string("\xc3\xa1"),      string("\xc3\xa9"),  string("\xc3\xad"),        string("\xc3\xb3"),      string("\xc3\xba"),    none,             none);
    res->add_row("empty",   string(),            string(),               string(),                string(),            string(),                  string(),                string(),              none,             string());
    // clang-format on
    return table_ptr(std::move(res));
}

// JSON
struct json_row
{
    std::string id;
    optional<std::string> field_json;
};
BOOST_DESCRIBE_STRUCT(json_row, (), (id, field_json))

static table_ptr types_json()
{
    // MariaDB doesn't have a dedicated column type, so there is a difference in metadata.
    // Values should be the same, though.
    auto res = make_table<json_row>("types_json");
    res->add_meta("field_json", get_server_features().json_type ? column_type::json : column_type::text);

    using std::string;

    // clang-format off
    res->add_row("regular",        string(R"([null, 42, false, "abc", {"key": "value"}])"));
    res->add_row("unicode_escape", string(R"(["\u0000value\u0000"])"));
    res->add_row("utf8",           string("[\"adi\xc3\xb3s\"]"));
    res->add_row("empty",          string("{}"));
    // clang-format on
    return table_ptr(std::move(res));
}

// Binary types
struct binary_row
{
    std::string id;
    optional<blob> field_binary;
    optional<blob> field_varbinary;
    optional<blob> field_tinyblob;
    optional<blob> field_blob;
    optional<blob> field_mediumblob;
    optional<blob> field_longblob;
};
BOOST_DESCRIBE_STRUCT(
    binary_row,
    (),
    (id, field_binary, field_varbinary, field_tinyblob, field_blob, field_mediumblob, field_longblob)
)

static table_ptr types_binary()
{
    auto res = make_table<binary_row>("types_binary");
    res->add_meta("field_binary", column_type::binary);
    res->add_meta("field_varbinary", column_type::varbinary);
    res->add_meta("field_tinyblob", column_type::blob);
    res->add_meta("field_blob", column_type::blob);
    res->add_meta("field_mediumblob", column_type::blob);
    res->add_meta("field_longblob", column_type::blob);

    // clang-format off
    res->add_row("regular",  makeb("\0_binary\0\0"),          makeb("\0_varbinary"),  makeb("\0_tinyblob"),  makeb("\0_blob"),  makeb("\0_mediumblob"), makeb("\0_longblob"));
    res->add_row("nonascii", makeb("\0\xff\0\0\0\0\0\0\0\0"), makeb("\1\xfe"),        makeb("\2\xfd"),       makeb("\3\xfc"),   makeb("\4\xfb"),        makeb("\5\xfa"));
    res->add_row("empty",    makeb("\0\0\0\0\0\0\0\0\0\0"),   blob(),                 blob(),                blob(),            blob(),                 blob());
    // clang-format on
    return table_ptr(std::move(res));
}

// These types don't have a better representation, and we represent
// them as strings or binary
struct not_implemented_row
{
    std::string id;
    optional<std::string> field_decimal;
    optional<blob> field_geometry;
};
BOOST_DESCRIBE_STRUCT(not_implemented_row, (), (id, field_decimal, field_geometry))

static table_ptr types_not_implemented()
{
    auto res = make_table<not_implemented_row>("types_not_implemented");
    res->add_meta("field_decimal", column_type::decimal);
    res->add_meta("field_geometry", column_type::geometry);

    blob geometry_value{0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                        0x00, 0x00, 0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40};

    res->add_row("regular", std::string("300"), geometry_value);
    return table_ptr(std::move(res));
}

// Tests for certain metadata flags
struct flags_row
{
    std::string id;
    optional<datetime> field_timestamp;
    int field_primary_key;
    std::string field_not_null;
    optional<int> field_unique;
    optional<int> field_indexed;
};
BOOST_DESCRIBE_STRUCT(
    flags_row,
    (),
    (id, field_timestamp, field_primary_key, field_not_null, field_unique, field_indexed)
)

static table_ptr types_flags()
{
    auto res = make_table<flags_row>("types_flags");
    res->add_meta(
        "field_timestamp",
        column_type::timestamp,
        flagsvec{&metadata::is_set_to_now_on_update},
        0,
        flags_unsigned
    );
    res->add_meta(
        "field_primary_key",
        column_type::int_,
        flagsvec{&metadata::is_primary_key, &metadata::is_not_null, &metadata::is_auto_increment}
    );
    res->add_meta("field_not_null", column_type::char_, flagsvec{&metadata::is_not_null});
    res->add_meta("field_unique", column_type::int_, flagsvec{&metadata::is_unique_key});
    res->add_meta("field_indexed", column_type::int_, flagsvec{&metadata::is_multiple_key});

    res->add_row("default", none, 50, "char", 21, 42);
    return table_ptr(std::move(res));
}

std::vector<table_ptr> make_all_tables()
{
    std::vector<table_ptr> res;
    res.push_back(types_tinyint());
    res.push_back(types_smallint());
    res.push_back(types_mediumint());
    res.push_back(types_int());
    res.push_back(types_bigint());
    res.push_back(types_year());
    res.push_back(types_bool());
    res.push_back(types_bit());
    res.push_back(types_float());
    res.push_back(types_double());
    res.push_back(types_date());
    res.push_back(types_datetime());
    res.push_back(types_timestamp());
    res.push_back(types_time());
    res.push_back(types_string());
    res.push_back(types_json());
    res.push_back(types_binary());
    res.push_back(types_not_implemented());
    res.push_back(types_flags());
    return res;
}

const std::vector<table_ptr>& all_tables()
{
    static std::vector<table_ptr> res = make_all_tables();
    return res;
}

BOOST_FIXTURE_TEST_CASE(query_read, database_types_fixture)
{
    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table->name)
        {
            // Execute the query
            results result;
            conn.async_execute(table->select_sql(), result, as_netresult).validate_no_error();

            // Validate the received contents
            validate_meta(result.meta(), table->metas);
            table->validate_rows(result.rows());
        }
    }
}

BOOST_FIXTURE_TEST_CASE(sql_format_query_write, database_types_fixture)
{
    start_transaction();

    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table->name)
        {
            // Remove all contents from the table
            results result;
            conn.async_execute(table->delete_sql(), result, as_netresult).validate_no_error();

            // Insert all the contents again
            conn.async_execute(table->insert_sql(), result, as_netresult).validate_no_error();

            // Query them again and verify the insertion was okay
            conn.async_execute(table->select_sql(), result, as_netresult).validate_no_error();
            validate_meta(result.meta(), table->metas);
            table->validate_rows(result.rows());
        }
    }
}

BOOST_FIXTURE_TEST_CASE(statement_read, database_types_fixture)
{
    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table->name)
        {
            // Prepare the statement
            auto stmt = conn.async_prepare_statement(table->select_sql(), as_netresult).get();

            // Execute it with the provided parameters
            results result;
            conn.async_execute(stmt.bind(), result, as_netresult).validate_no_error();

            // Validate the received contents
            validate_meta(result.meta(), table->metas);
            table->validate_rows(result.rows());
        }
    }
}

BOOST_FIXTURE_TEST_CASE(statement_write, database_types_fixture)
{
    start_transaction();

    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table->name)
        {
            // Prepare the statements
            auto insert_stmt = conn.async_prepare_statement(table->insert_sql_stmt(), as_netresult).get();
            auto query_stmt = conn.async_prepare_statement(table->select_sql(), as_netresult).get();

            // Remove all contents from the table
            results result;
            conn.async_execute(table->delete_sql(), result, as_netresult).validate_no_error();

            // Insert all the contents again
            for (const auto& row : table->rws)
            {
                conn.async_execute(insert_stmt.bind(row.begin(), row.end()), result, as_netresult)
                    .validate_no_error();
            }

            // Query them again and verify the insertion was okay
            conn.async_execute(query_stmt.bind(), result, as_netresult).validate_no_error();
            validate_meta(result.meta(), table->metas);
            table->validate_rows(result.rows());
        }
    }
}

#ifdef BOOST_MYSQL_CXX14
BOOST_FIXTURE_TEST_CASE(static_interface, database_types_fixture)
{
    for (const auto& table : all_tables())
    {
        BOOST_TEST_CONTEXT(table->name)
        {
            // All the test is type-specific
            table->select_static(conn);
        }
    }
}
#endif

BOOST_AUTO_TEST_SUITE_END()  // test_database_types
