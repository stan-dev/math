//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/datetime.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "test_common/stringize.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_datetime)

// Helpers
datetime from_timestamp(std::int64_t micros_since_epoch)
{
    return datetime(datetime::time_point(datetime::time_point::duration(micros_since_epoch)));
}

#ifdef BOOST_MYSQL_HAS_LOCAL_TIME
datetime from_local_timestamp(std::int64_t micros_since_epoch)
{
    return datetime(datetime::local_time_point(datetime::local_time_point::duration(micros_since_epoch)));
}
#endif

BOOST_AUTO_TEST_CASE(default_ctor)
{
    datetime d;
    BOOST_TEST(d.year() == 0u);
    BOOST_TEST(d.month() == 0u);
    BOOST_TEST(d.day() == 0u);
    BOOST_TEST(d.hour() == 0u);
    BOOST_TEST(d.minute() == 0u);
    BOOST_TEST(d.second() == 0u);
    BOOST_TEST(d.microsecond() == 0u);
    BOOST_TEST(!d.valid());
}

// Also verifies that all the passed datetimes are considered valid
BOOST_AUTO_TEST_CASE(from_to_time_point)
{
    struct
    {
        const char* name;
        std::int64_t micros_since_epoch;
        datetime d;
    } test_cases[] = {
        {"date",                   1272758400000000,   datetime(2010, 5,  2,  0,  0,  0,  0)     },
        {"date_leap4",             1078012800000000,   datetime(2004, 2,  29, 0,  0,  0,  0)     },
        {"date_leap400",           951782400000000,    datetime(2000, 2,  29, 0,  0,  0,  0)     },
        {"u",                      1272758400123456,   datetime(2010, 5,  2,  0,  0,  0,  123456)},
        {"s",                      1272758450000000,   datetime(2010, 5,  2,  0,  0,  50, 0)     },
        {"m",                      1272758460000000,   datetime(2010, 5,  2,  0,  1,  0,  0)     },
        {"hs",                     1272841250000000,   datetime(2010, 5,  2,  23, 0,  50, 0)     },
        {"ms",                     1272758510000000,   datetime(2010, 5,  2,  0,  1,  50, 0)     },
        {"hu",                     1272841200123456,   datetime(2010, 5,  2,  23, 0,  0,  123456)},
        {"mu",                     1272758460123456,   datetime(2010, 5,  2,  0,  1,  0,  123456)},
        {"hmu",                    1272841260123456,   datetime(2010, 5,  2,  23, 1,  0,  123456)},
        {"su",                     1272758450123456,   datetime(2010, 5,  2,  0,  0,  50, 123456)},
        {"hsu",                    1272841250123456,   datetime(2010, 5,  2,  23, 0,  50, 123456)},
        {"msu",                    1272758510123456,   datetime(2010, 5,  2,  0,  1,  50, 123456)},
        {"h",                      1272841200000000,   datetime(2010, 5,  2,  23, 0,  0,  0)     },
        {"hm",                     1272841260000000,   datetime(2010, 5,  2,  23, 1,  0,  0)     },
        {"hms",                    1272841310000000,   datetime(2010, 5,  2,  23, 1,  50, 0)     },
        {"hmsu",                   1272841310123456,   datetime(2010, 5,  2,  23, 1,  50, 123456)},
        {"hmsu_minfraction",       1262307661000001,   datetime(2010, 1,  1,  1,  1,  1,  1)     },
        {"hmsu_halffraction1",     1276601369499999,   datetime(2010, 6,  15, 11, 29, 29, 499999)},
        {"hmsu_halffraction2",     1276605030500000,   datetime(2010, 6,  15, 12, 30, 30, 500000)},
        {"hmsu_maxfraction2",      1293839999999999,   datetime(2010, 12, 31, 23, 59, 59, 999999)},

        {"neg_date",               -86400000000,       datetime(1969, 12, 31, 0,  0,  0,  0)     },
        {"neg_date_2",             -3 * 86400000000,   datetime(1969, 12, 29, 0,  0,  0,  0)     },
        {"neg_date_leap4",         -11544768000000000, datetime(1604, 2,  29, 0,  0,  0,  0)     },
        {"neg_date_leap400",       -11670998400000000, datetime(1600, 2,  29, 0,  0,  0,  0)     },
        {"neg_u",                  -84239999876544,    datetime(1967, 5,  2,  0,  0,  0,  123456)},
        {"neg_s",                  -84239950000000,    datetime(1967, 5,  2,  0,  0,  50, 0)     },
        {"neg_m",                  -84239940000000,    datetime(1967, 5,  2,  0,  1,  0,  0)     },
        {"neg_hs",                 -84157150000000,    datetime(1967, 5,  2,  23, 0,  50, 0)     },
        {"neg_ms",                 -84239890000000,    datetime(1967, 5,  2,  0,  1,  50, 0)     },
        {"neg_hu",                 -84157199876544,    datetime(1967, 5,  2,  23, 0,  0,  123456)},
        {"neg_mu",                 -84239939876544,    datetime(1967, 5,  2,  0,  1,  0,  123456)},
        {"neg_hmu",                -84157139876544,    datetime(1967, 5,  2,  23, 1,  0,  123456)},
        {"neg_su",                 -84239949876544,    datetime(1967, 5,  2,  0,  0,  50, 123456)},
        {"neg_hsu",                -84157149876544,    datetime(1967, 5,  2,  23, 0,  50, 123456)},
        {"neg_msu",                -84239889876544,    datetime(1967, 5,  2,  0,  1,  50, 123456)},
        {"neg_h",                  -84157200000000,    datetime(1967, 5,  2,  23, 0,  0,  0)     },
        {"neg_hm",                 -84157140000000,    datetime(1967, 5,  2,  23, 1,  0,  0)     },
        {"neg_hms",                -84157090000000,    datetime(1967, 5,  2,  23, 1,  50, 0)     },
        {"neg_hmsu",               -84157089876544,    datetime(1967, 5,  2,  23, 1,  50, 123456)},
        {"neg_hmsu_minfraction",   -1893452338999999,  datetime(1910, 1,  1,  1,  1,  1,  1)     },
        {"neg_hmsu_halffraction1", -3141376230500001,  datetime(1870, 6,  15, 11, 29, 29, 499999)},
        {"neg_hmsu_halffraction2", -5697430169500000,  datetime(1789, 6,  15, 12, 30, 30, 500000)},
        {"neg_hmsu_maxfraction2",  -47114438400000001, datetime(476,  12, 31, 23, 59, 59, 999999)},

        {"epoch",                  0,                  datetime(1970, 1,  1,  0,  0,  0,  0)     },
        {"min",                    -62167219200000000, datetime(0,    1,  1,  0,  0,  0,  0)     },
        {"max",                    253402300799999999, datetime(9999, 12, 31, 23, 59, 59, 999999)},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // valid
            BOOST_TEST(tc.d.valid());

            // From time_point
            datetime actual_dt = from_timestamp(tc.micros_since_epoch);
            BOOST_TEST(actual_dt == tc.d);

            // To time_point
            auto tp = tc.d.get_time_point();
            BOOST_TEST(tp.time_since_epoch().count() == tc.micros_since_epoch);
        }
    }
}

// All cases that are accepeted by datetime but don't represent actual datetimes
BOOST_AUTO_TEST_CASE(valid_false)
{
    struct
    {
        const char* name;
        datetime d;
    } test_cases[] = {
        {"date_yregular_invalid_date",             datetime(2020,   11,   31,   0,    0,    0,    0)         },
        {"date_yregular_invalid_date_leap100",     datetime(1900,   2,    29,   0,    0,    0,    0)         },
        {"date_yregular_invalid_date_leapregular", datetime(1999,   2,    29,   0,    0,    0,    0)         },
        {"date_yregular_mregular_dzero",           datetime(2020,   10,   0,    0,    0,    0,    0)         },
        {"date_yregular_mzero_dregular",           datetime(2020,   0,    10,   0,    0,    0,    0)         },
        {"date_yregular_mzero_dzero",              datetime(2020,   0,    0,    0,    0,    0,    0)         },
        {"date_yzero_invalid_date",                datetime(0,      11,   31,   0,    0,    0,    0)         },
        {"date_yzero_mregular_dzero",              datetime(0,      10,   0,    0,    0,    0,    0)         },
        {"date_yzero_mzero_dregular",              datetime(0,      0,    10,   0,    0,    0,    0)         },
        {"date_zero",                              datetime(0,      0,    0,    0,    0,    0,    0)         },
        {"hms_yregular_invalid_date",              datetime(2020,   11,   31,   10,   20,   30,   0)         },
        {"hms_yregular_invalid_date_leap100",      datetime(1900,   2,    29,   10,   20,   30,   0)         },
        {"hms_yregular_invalid_date_leapregular",  datetime(1999,   2,    29,   10,   20,   30,   0)         },
        {"hms_yregular_mregular_dzero",            datetime(2020,   10,   0,    10,   20,   30,   0)         },
        {"hms_yregular_mzero_dregular",            datetime(2020,   0,    10,   10,   20,   30,   0)         },
        {"hms_yregular_mzero_dzero",               datetime(2020,   0,    0,    10,   20,   30,   0)         },
        {"hms_yzero_invalid_date",                 datetime(0,      11,   31,   10,   20,   30,   0)         },
        {"hms_yzero_mregular_dzero",               datetime(0,      10,   0,    10,   20,   30,   0)         },
        {"hms_yzero_mzero_dregular",               datetime(0,      0,    10,   10,   20,   30,   0)         },
        {"hms_zero",                               datetime(0,      0,    0,    10,   20,   30,   0)         },
        {"hmsu_yregular_invalid_date",             datetime(2020,   11,   31,   10,   20,   30,   999999)    },
        {"hmsu_yregular_invalid_date_leap100",     datetime(1900,   2,    29,   10,   20,   30,   999999)    },
        {"hmsu_yregular_invalid_date_leapregular", datetime(1999,   2,    29,   10,   20,   30,   999999)    },
        {"hmsu_yregular_mregular_dzero",           datetime(2020,   10,   0,    10,   20,   30,   999999)    },
        {"hmsu_yregular_mzero_dregular",           datetime(2020,   0,    10,   10,   20,   30,   999999)    },
        {"hmsu_yregular_mzero_dzero",              datetime(2020,   0,    0,    10,   20,   30,   999999)    },
        {"hmsu_yzero_invalid_date",                datetime(0,      11,   31,   10,   20,   30,   999999)    },
        {"hmsu_yzero_mregular_dzero",              datetime(0,      10,   0,    10,   20,   30,   999999)    },
        {"hmsu_yzero_mzero_dregular",              datetime(0,      0,    10,   10,   20,   30,   999999)    },
        {"hmsu_zero",                              datetime(0,      0,    0,    10,   20,   30,   999999)    },
        {"invalid_hour",                           datetime(2010,   5,    22,   24,   0,    0,    0)         },
        {"invalid_minute",                         datetime(2010,   5,    22,   23,   60,   0,    0)         },
        {"invalid_second",                         datetime(2010,   5,    22,   23,   59,   60,   0)         },
        {"invalid_microsecond",                    datetime(2010,   5,    22,   23,   59,   59,   1000000)   },
        {"max_year",                               datetime(0xffff, 5,    2,    23,   1,    50,   123456)    },
        {"max_month",                              datetime(2010,   0xff, 2,    23,   1,    50,   123456)    },
        {"max_day",                                datetime(2010,   5,    0xff, 23,   1,    50,   123456)    },
        {"max_hour",                               datetime(2010,   5,    2,    0xff, 1,    50,   123456)    },
        {"max_minute",                             datetime(2010,   5,    2,    23,   0xff, 50,   123456)    },
        {"max_second",                             datetime(2010,   5,    2,    23,   1,    0xff, 123456)    },
        {"max_microsecond",                        datetime(2010,   5,    2,    23,   1,    50,   0xffffffff)},
        {"max",                                    datetime(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff)},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name) { BOOST_TEST(!tc.d.valid()); }
    }
}

BOOST_AUTO_TEST_CASE(ctor_from_time_point_invalid)
{
    BOOST_CHECK_THROW(from_timestamp(253402300799999999 + 1), std::out_of_range);
    BOOST_CHECK_THROW(from_timestamp(-62167219200000000 - 1), std::out_of_range);
    BOOST_CHECK_THROW(from_timestamp((std::numeric_limits<std::int64_t>::max)()), std::out_of_range);
    BOOST_CHECK_THROW(from_timestamp((std::numeric_limits<std::int64_t>::min)()), std::out_of_range);
}

#ifdef BOOST_MYSQL_HAS_LOCAL_TIME
BOOST_AUTO_TEST_CASE(ctor_from_local_time_point)
{
    BOOST_TEST(from_local_timestamp(253402300799999999) == datetime(9999, 12, 31, 23, 59, 59, 999999));
    BOOST_TEST(from_local_timestamp(1715966176806454) == datetime(2024, 5, 17, 17, 16, 16, 806454));
    BOOST_TEST(from_local_timestamp(-62167219200000000) == datetime(0, 1, 1));
}

BOOST_AUTO_TEST_CASE(ctor_from_local_time_point_invalid)
{
    BOOST_CHECK_THROW(from_local_timestamp(253402300799999999 + 1), std::out_of_range);
    BOOST_CHECK_THROW(from_local_timestamp(-62167219200000000 - 1), std::out_of_range);
    BOOST_CHECK_THROW(from_local_timestamp((std::numeric_limits<std::int64_t>::max)()), std::out_of_range);
    BOOST_CHECK_THROW(from_local_timestamp((std::numeric_limits<std::int64_t>::min)()), std::out_of_range);
}
#endif

// spotcheck, uses the same routines as get_time_point
BOOST_AUTO_TEST_CASE(as_time_point)
{
    datetime d1(2010, 12, 31, 23, 59, 59, 999999);
    BOOST_TEST(d1.as_time_point().time_since_epoch().count() == 1293839999999999);

    datetime d2(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff);
    BOOST_CHECK_THROW(d2.as_time_point(), std::invalid_argument);
}

#ifdef BOOST_MYSQL_HAS_LOCAL_TIME
// spotcheck, uses the same routines as get_time_point
BOOST_AUTO_TEST_CASE(as_local_time_point)
{
    datetime d1(2010, 12, 31, 23, 59, 59, 999999);
    BOOST_TEST(d1.as_local_time_point().time_since_epoch().count() == 1293839999999999);

    datetime d2(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff);
    BOOST_CHECK_THROW(d2.as_local_time_point(), std::invalid_argument);
}

// spotcheck, uses the same routines as get_time_point
BOOST_AUTO_TEST_CASE(get_local_time_point)
{
    datetime d(2010, 12, 31, 23, 59, 59, 999999);
    BOOST_TEST(d.get_local_time_point().time_since_epoch().count() == 1293839999999999);
}
#endif

BOOST_AUTO_TEST_CASE(operator_equals)
{
    datetime maxdt(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff);

    // clang-format off
    struct
    {
        const char* name;
        datetime d1;
        datetime d2;
        bool equal;
    } test_cases[] = {
        {"equal", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2020, 2, 29, 1, 2, 10, 9), true},
        {"equal_invalid", datetime(0, 0, 0, 0, 0, 0, 0), datetime(0, 0, 0, 0, 0, 0, 0), true},
        {"equal_max", maxdt, maxdt, true},
        {"ne_year", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2019, 2, 29, 1, 2, 10, 9), false},
        {"ne_month", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2020, 1, 29, 1, 2, 10, 9), false},
        {"ne_day", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2020, 2, 28, 1, 2, 10, 9), false},
        {"ne_hour", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2020, 2, 29, 2, 2, 10, 9), false},
        {"ne_min", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2020, 2, 29, 1, 3, 10, 9), false},
        {"ne_sec", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2020, 2, 29, 1, 2, 11, 9), false},
        {"ne_micro", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2020, 2, 29, 1, 2, 10, 98), false},
        {"ne_all", datetime(2020, 2, 29, 1, 2, 10, 9), datetime(2019, 5, 10, 0xff, 21, 101, 98), false},
    };
    // clang-format on

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST((tc.d1 == tc.d2) == tc.equal);
            BOOST_TEST((tc.d2 == tc.d1) == tc.equal);
            BOOST_TEST((tc.d1 != tc.d2) == !tc.equal);
            BOOST_TEST((tc.d2 != tc.d1) == !tc.equal);
        }
    }
}

// operator<< is implemented in terms of datetime_to_string (see dt_to_string.hpp)
BOOST_AUTO_TEST_CASE(operator_stream)
{
    BOOST_TEST(stringize(datetime(2023, 1, 2, 12, 10, 1, 0)) == "2023-01-02 12:10:01.000000");
    BOOST_TEST(stringize(datetime(2022, 12, 31)) == "2022-12-31 00:00:00.000000");
    BOOST_TEST(stringize(datetime(2020, 3, 2, 23, 59, 59, 12345)) == "2020-03-02 23:59:59.012345");
    BOOST_TEST(stringize(datetime()) == "0000-00-00 00:00:00.000000");
    BOOST_TEST(
        stringize(datetime(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff)) ==
        "65535-255-255 255:255:255.4294967295"
    );
}

BOOST_AUTO_TEST_CASE(now)
{
    auto d = datetime::now();
    BOOST_TEST_REQUIRE(d.valid());
    BOOST_TEST(d.year() > 2020u);
    BOOST_TEST(d.year() < 2100u);
}

// Make sure constxpr can actually be used in a constexpr context
BOOST_AUTO_TEST_CASE(constexpr_fns_cxx11)
{
    constexpr datetime d0{};
    static_assert(!d0.valid(), "");
    static_assert(d0.year() == 0, "");
    static_assert(d0.month() == 0, "");
    static_assert(d0.day() == 0, "");
    static_assert(d0.hour() == 0, "");
    static_assert(d0.minute() == 0, "");
    static_assert(d0.second() == 0, "");
    static_assert(d0.microsecond() == 0, "");

    constexpr datetime d1(2020, 10, 1, 23, 20, 59, 123456);
    static_assert(d1.valid(), "");
    static_assert(d1.year() == 2020, "");
    static_assert(d1.month() == 10, "");
    static_assert(d1.day() == 1, "");
    static_assert(d1.hour() == 23, "");
    static_assert(d1.minute() == 20, "");
    static_assert(d1.second() == 59, "");
    static_assert(d1.microsecond() == 123456, "");

    static_assert(d0 == datetime(), "");
    static_assert(d0 != d1, "");
}

#ifndef BOOST_NO_CXX14_CONSTEXPR
BOOST_AUTO_TEST_CASE(constexpr_fns_cxx14)
{
    constexpr datetime d0(datetime::time_point(datetime::time_point::duration(1272841310123456)));
    static_assert(d0 == datetime(2010, 5, 2, 23, 1, 50, 123456), "");

    constexpr auto tp1 = d0.get_time_point();
    static_assert(tp1.time_since_epoch().count() == 1272841310123456, "");

    constexpr auto tp2 = d0.as_time_point();
    static_assert(tp2 == tp1, "");
}
#endif

BOOST_AUTO_TEST_SUITE_END()
