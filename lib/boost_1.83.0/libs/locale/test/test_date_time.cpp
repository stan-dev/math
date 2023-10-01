//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/date_time.hpp>
#include <boost/locale/formatting.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/localization_backend.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <limits>

#ifdef BOOST_LOCALE_WITH_ICU
#    include <unicode/uversion.h>
#    define BOOST_LOCALE_ICU_VERSION (U_ICU_VERSION_MAJOR_NUM * 100 + U_ICU_VERSION_MINOR_NUM)
#else
#    define BOOST_LOCALE_ICU_VERSION 0
#endif

#ifdef BOOST_MSVC
#    pragma warning(disable : 4244) // loose data
#endif

#define TEST_EQ_FMT(t, X)    \
    empty_stream(ss) << (t); \
    test_eq_impl(ss.str(), X, #t "==" #X, __LINE__)

// Very simple container for a part of the tests. Counts its instances
struct mock_calendar : public boost::locale::abstract_calendar {
    using period_mark = boost::locale::period::marks::period_mark;

    mock_calendar() : time(0) { ++num_instances; }
    mock_calendar(const mock_calendar& other) : time(other.time) { ++num_instances; }
    ~mock_calendar() { --num_instances; }

    abstract_calendar* clone() const override { return new mock_calendar(*this); }
    void set_value(period_mark, int) override {}                        // LCOV_EXCL_LINE
    void normalize() override {}                                        // LCOV_EXCL_LINE
    int get_value(period_mark, value_type) const override { return 0; } // LCOV_EXCL_LINE
    void set_time(const boost::locale::posix_time&) override {}         // LCOV_EXCL_LINE
    boost::locale::posix_time get_time() const override { return {}; }  // LCOV_EXCL_LINE
    double get_time_ms() const override { return time; }
    void set_option(calendar_option_type, int) override {}                             // LCOV_EXCL_LINE
    int get_option(calendar_option_type) const override { return 0; }                  // LCOV_EXCL_LINE
    void adjust_value(period_mark, update_type, int) override {}                       // LCOV_EXCL_LINE
    int difference(const abstract_calendar&, period_mark) const override { return 0; } // LCOV_EXCL_LINE
    void set_timezone(const std::string&) override {}
    std::string get_timezone() const override { return "mock TZ"; }
    bool same(const abstract_calendar* other) const override
    {
        return dynamic_cast<const mock_calendar*>(other) != nullptr;
    }

    static int num_instances;
    double time;
};
int mock_calendar::num_instances = 0;
struct mock_calendar_facet : boost::locale::calendar_facet {
    boost::locale::abstract_calendar* create_calendar() const override { return proto_cal.clone(); }
    mock_calendar proto_cal;
};

void test_main(int /*argc*/, char** /*argv*/)
{
    using namespace boost::locale;
    using namespace boost::locale::period;
    std::unique_ptr<calendar> mock_cal;
    {
        auto* cal_facet = new mock_calendar_facet;
        std::locale old_loc = std::locale::global(std::locale(std::locale(), cal_facet));
        mock_calendar::num_instances = 0;
        {
            cal_facet->proto_cal.time = 42 * 1e3;
            date_time t1;
            TEST_EQ(t1.time(), 42);
            TEST_EQ(t1.timezone(), "mock TZ");
            TEST_EQ(mock_calendar::num_instances, 1);
            cal_facet->proto_cal.time = 99 * 1e3;
            date_time t2;
            TEST_EQ(t2.time(), 99);
            TEST_EQ(mock_calendar::num_instances, 2);
            // Copy construct
            date_time t3 = t1;
            TEST_EQ(t1.time(), 42);
            TEST_EQ(t2.time(), 99);
            TEST_EQ(t3.time(), 42);
            TEST_EQ(mock_calendar::num_instances, 3);
            // Copy assign
            t3 = t2;
            TEST_EQ(t3.time(), 99);
            TEST_EQ(mock_calendar::num_instances, 3); // No new
            {
                // Move construct
                date_time t4 = std::move(t1);
                TEST_EQ(t4.time(), 42);
                TEST_EQ(mock_calendar::num_instances, 3); // No new
                // Move assign
                t2 = std::move(t4);
                TEST_EQ(t2.time(), 42);
                TEST_LE(mock_calendar::num_instances, 3); // maybe destroy old t2
            }
            // Unchanged after t4 (or old t2) is destroyed
            TEST_EQ(t2.time(), 42);
            TEST_EQ(mock_calendar::num_instances, 2);
            // Self move, via reference to avoid triggering a compiler warning
            date_time& t2_ref = t2;
            t2_ref = std::move(t2);
            TEST_EQ(t2.time(), 42);
            TEST_EQ(mock_calendar::num_instances, 2);
        }
        TEST_EQ(mock_calendar::num_instances, 0); // No leaks
        mock_cal.reset(new calendar());
        std::locale::global(old_loc);
    }
    for(const std::string& backend_name : boost::locale::localization_backend_manager::global().get_all_backends()) {
        std::cout << "Testing for backend: " << backend_name << std::endl;
        boost::locale::localization_backend_manager tmp_backend = boost::locale::localization_backend_manager::global();
        tmp_backend.select(backend_name);
        boost::locale::localization_backend_manager::global(tmp_backend);

        boost::locale::generator g;
        std::locale loc = g("en_US.UTF-8");
        {
            using boost::locale::abstract_calendar;
            std::unique_ptr<abstract_calendar> cal(
              std::use_facet<boost::locale::calendar_facet>(loc).create_calendar());
            TEST_THROWS(cal->set_option(abstract_calendar::is_gregorian, 0), boost::locale::date_time_error);
            TEST_THROWS(cal->set_option(abstract_calendar::is_dst, 0), boost::locale::date_time_error);
        }

        {
            std::locale::global(loc);

            const std::string tz = "GMT";
            time_zone::global(tz);
            // A call returns the old tz
            TEST_EQ(time_zone::global("GMT+01:00"), tz);
            TEST_EQ(time_zone::global(tz), "GMT+01:00");
            calendar cal(loc, tz);
            TEST(cal.get_locale() == loc);
            TEST_EQ(cal.get_time_zone(), tz);

            TEST(calendar() == cal);
            TEST(calendar(loc) == cal);
            TEST(calendar(tz) == cal);
            {
                const std::string tz2 = "GMT+01:00";
                const std::locale loc2 = g("ru_RU.UTF-8");
                const calendar cal_tz2(loc, "GMT+01:00");
                const calendar cal_loc2(loc2);
                TEST(cal_tz2 != cal);
                TEST(cal_loc2 != cal);
                calendar cal_tmp(cal);
                TEST(cal_tmp == cal);
                TEST(cal_tmp != cal_tz2);
                cal_tmp = cal_tz2;
                TEST(cal_tmp == cal_tz2);
                TEST_EQ(cal_tmp.get_time_zone(), tz2);
                TEST(cal_tmp.get_locale() == loc);
                TEST(cal_tmp != cal_loc2);
                cal_tmp = cal_loc2;
                TEST(cal_tmp == cal_loc2);
                TEST_EQ(cal_tmp.get_time_zone(), tz);
                TEST(cal_tmp.get_locale() == loc2);
            }
            {
                calendar cal2;
                TEST(cal2 != *mock_cal);
                cal2 = *mock_cal;
                TEST(cal2 == *mock_cal);
            }

            TEST_EQ(cal.minimum(month()), 0);
            TEST_EQ(cal.maximum(month()), 11);
            TEST_EQ(cal.minimum(day()), 1);
            TEST_EQ(cal.greatest_minimum(day()), 1);
            TEST_EQ(cal.least_maximum(day()), 28);
            TEST_EQ(cal.maximum(day()), 31);
            TEST(cal.is_gregorian());

            TEST_EQ(calendar(g("ar_EG.UTF-8")).first_day_of_week(), 7);
            TEST_EQ(calendar(g("he_IL.UTF-8")).first_day_of_week(), 1);
            TEST_EQ(calendar(g("ru_RU.UTF-8")).first_day_of_week(), 2);

            std::ostringstream ss;
            ss.imbue(loc);
            ss << boost::locale::as::time_zone(tz);

            const time_t one_h = 60 * 60;
            const time_t a_date = 24 * one_h * (31 + 4);     // Feb 5th
            const time_t a_time = 15 * one_h + 60 * 33 + 13; // 15:33:13
            const time_t a_datetime = a_date + a_time;

            const date_time tp_5_feb_1970_153313 = date_time(a_datetime); // 5th Feb 1970 15:33:13
            TEST_EQ(tp_5_feb_1970_153313.timezone(), tz);
            ss << as::ftime("%Y-%m-%d");
            TEST_EQ_FMT(tp_5_feb_1970_153313, "1970-02-05");
            ss << as::ftime("%Y-%m-%d %H:%M:%S");
            TEST_EQ_FMT(tp_5_feb_1970_153313, "1970-02-05 15:33:13");

            // Test set()
            date_time time_point = tp_5_feb_1970_153313;
            time_point.set(year(), 1990);
            TEST_EQ_FMT(time_point, "1990-02-05 15:33:13");
            time_point.set(month(), 5);
            TEST_EQ_FMT(time_point, "1990-06-05 15:33:13");
            time_point.set(day(), 9);
            TEST_EQ_FMT(time_point, "1990-06-09 15:33:13");
            time_point.set(hour(), 11);
            TEST_EQ_FMT(time_point, "1990-06-09 11:33:13");
            time_point.set(minute(), 42);
            TEST_EQ_FMT(time_point, "1990-06-09 11:42:13");
            time_point.set(second(), 24);
            TEST_EQ_FMT(time_point, "1990-06-09 11:42:24");
            time_point.set(am_pm(), 1);
            TEST_EQ_FMT(time_point, "1990-06-09 23:42:24");
            // Overflow day of month
            time_point.set(day(), time_point.maximum(day()) + 1);
            TEST_EQ_FMT(time_point, "1990-07-01 23:42:24");

            // Same via assignment
            time_point = tp_5_feb_1970_153313;
            time_point = year(1990);
            TEST_EQ_FMT(time_point, "1990-02-05 15:33:13");
            time_point = month(5);
            TEST_EQ_FMT(time_point, "1990-06-05 15:33:13");
            time_point = day(9);
            TEST_EQ_FMT(time_point, "1990-06-09 15:33:13");
            time_point = hour(11);
            TEST_EQ_FMT(time_point, "1990-06-09 11:33:13");
            time_point = minute(42);
            TEST_EQ_FMT(time_point, "1990-06-09 11:42:13");
            time_point = second(24);
            TEST_EQ_FMT(time_point, "1990-06-09 11:42:24");
            time_point = am_pm(1);
            TEST_EQ_FMT(time_point, "1990-06-09 23:42:24");
            // Overflow day of month
            time_point = day(time_point.maximum(day()) + 1);
            TEST_EQ_FMT(time_point, "1990-07-01 23:42:24");
            // All at once
            time_point = year(1989) + month(2) + day(5) + hour(7) + minute(9) + second(11);
            TEST_EQ_FMT(time_point, "1989-03-05 07:09:11");
            // Partials:
            time_point = year(1970) + february() + day(5);
            TEST_EQ_FMT(time_point, "1970-02-05 07:09:11");
            time_point = 3 * hour_12() + 1 * am_pm() + 33 * minute() + 13 * second();
            TEST_EQ_FMT(time_point, "1970-02-05 15:33:13");

            time_point = tp_5_feb_1970_153313;
            time_point += hour();
            TEST_EQ_FMT(time_point, "1970-02-05 16:33:13");

            TEST_EQ(time_point.minimum(day()), 1);
            TEST_EQ(time_point.maximum(day()), 28);

            time_point = tp_5_feb_1970_153313;
            time_point += year() * 2 + 1 * month();
            TEST_EQ_FMT(time_point, "1972-03-05 15:33:13");

            time_point = tp_5_feb_1970_153313;
            time_point -= minute();
            TEST_EQ_FMT(time_point, "1970-02-05 15:32:13");

            time_point = tp_5_feb_1970_153313;
            time_point <<= minute() * 30;
            TEST_EQ_FMT(time_point, "1970-02-05 15:03:13");

            time_point = tp_5_feb_1970_153313;
            time_point >>= minute(40);
            TEST_EQ_FMT(time_point, "1970-02-05 15:53:13");

            time_point = tp_5_feb_1970_153313;
            TEST_EQ((time_point + month()) / month(), 2);
            TEST_EQ(month(time_point + month(1)), 2);
            TEST_EQ(time_point / month(), 1);
            TEST_EQ((time_point - month()) / month(), 0);
            TEST_EQ(time_point / month(), 1);
            TEST_EQ((time_point << month()) / month(), 2);
            TEST_EQ(time_point / month(), 1);
            TEST_EQ((time_point >> month()) / month(), 0);
            TEST_EQ(time_point / month(), 1);

            // To subtract from the year, don't use 1970 which may be the lowest possible year
            const date_time tp_5_april_1990_153313 = (date_time(tp_5_feb_1970_153313) = (year(1990) + april()));
            TEST_EQ_FMT(tp_5_april_1990_153313, "1990-04-05 15:33:13");
            // Test each period
            TEST_EQ_FMT(tp_5_april_1990_153313 + year(2), "1992-04-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << year(2), "1992-04-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - year(10), "1980-04-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> year(10), "1980-04-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + month(2), "1990-06-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << month(2), "1990-06-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - month(1), "1990-03-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> month(1), "1990-03-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + day(2), "1990-04-07 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << day(2), "1990-04-07 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - day(3), "1990-04-02 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> day(3), "1990-04-02 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + hour(2), "1990-04-05 17:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << hour(2), "1990-04-05 17:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - hour(3), "1990-04-05 12:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> hour(3), "1990-04-05 12:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + minute(2), "1990-04-05 15:35:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << minute(2), "1990-04-05 15:35:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - minute(3), "1990-04-05 15:30:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> minute(3), "1990-04-05 15:30:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + second(2), "1990-04-05 15:33:15");
            TEST_EQ_FMT(tp_5_april_1990_153313 << second(2), "1990-04-05 15:33:15");
            TEST_EQ_FMT(tp_5_april_1990_153313 - second(2), "1990-04-05 15:33:11");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> second(2), "1990-04-05 15:33:11");
            // Difference between add and roll: The latter only changes the given field
            // So this tests what happens when going over/under the bound for each field
            TEST_EQ_FMT(tp_5_april_1990_153313 + month(12 + 2), "1991-06-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << month(12 + 2), "1990-06-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - month(12 * 3 + 1), "1987-03-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> month(12 * 3 + 1), "1990-03-05 15:33:13");
            // Check that possible int overflows get handled
            constexpr int max_full_years_in_months = (std::numeric_limits<int>::max() / 12) * 12;
            TEST_EQ_FMT(tp_5_april_1990_153313 >> month(max_full_years_in_months), "1990-04-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << month(max_full_years_in_months), "1990-04-05 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + day(30 + 2), "1990-05-07 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << day(30 + 2), "1990-04-07 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - day(10), "1990-03-26 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> day(10), "1990-04-25 15:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + hour(24 * 3 + 2), "1990-04-08 17:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << hour(24 * 3 + 2), "1990-04-05 17:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - hour(24 * 5 + 3), "1990-03-31 12:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> hour(24 * 5 + 3), "1990-04-05 12:33:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + minute(60 * 5 + 3), "1990-04-05 20:36:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 << minute(60 * 5 + 3), "1990-04-05 15:36:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 - minute(60 * 5 + 3), "1990-04-05 10:30:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> minute(60 * 5 + 3), "1990-04-05 15:30:13");
            TEST_EQ_FMT(tp_5_april_1990_153313 + second(60 * 3 + 2), "1990-04-05 15:36:15");
            TEST_EQ_FMT(tp_5_april_1990_153313 << second(60 * 3 + 2), "1990-04-05 15:33:15");
            TEST_EQ_FMT(tp_5_april_1990_153313 - second(60 * 5 + 2), "1990-04-05 15:28:11");
            TEST_EQ_FMT(tp_5_april_1990_153313 >> second(60 * 5 + 2), "1990-04-05 15:33:11");

            // Add a set of periods
            TEST_EQ_FMT(tp_5_feb_1970_153313 << (year(2) + month(3) - day(1) + hour(5) + minute(7) + second(9)),
                        "1972-05-04 20:40:22");
            TEST_EQ_FMT(tp_5_feb_1970_153313 + (year(2) + month(3) - day(1) + hour(5) + minute(7) + second(9)),
                        "1972-05-04 20:40:22");
            // std calendar can't go below 1970
            time_point = tp_5_feb_1970_153313;
            time_point = year(1972) + july();
            TEST_EQ_FMT(time_point, "1972-07-05 15:33:13");
            TEST_EQ_FMT(time_point >> (year(2) + month(3) - day(11) + hour(5) + minute(7) + second(9)),
                        "1970-04-16 10:26:04");
            TEST_EQ_FMT(time_point - (year(2) + month(3) - day(11) + hour(5) + minute(7) + second(9)),
                        "1970-04-16 10:26:04");

            time_point = tp_5_feb_1970_153313;
            TEST(time_point == tp_5_feb_1970_153313);
            TEST(!(time_point != tp_5_feb_1970_153313));
            TEST_EQ(time_point.get(hour()), 15);
            TEST_EQ(time_point / hour(), 15);
            TEST(time_point + year() != time_point);
            TEST(time_point - minute() <= time_point);
            TEST(time_point <= time_point);
            TEST(time_point + minute() >= time_point);
            TEST(time_point >= time_point);

            TEST(time_point < time_point + second());
            TEST(!(time_point < time_point - second()));
            TEST(time_point > time_point - second());
            TEST(!(time_point > time_point + second()));

            TEST_EQ(time_point.get(day()), 5);
            TEST_EQ(time_point.get(year()), 1970);

            TEST_EQ(time_point.get(era()), 1);
            TEST_EQ(time_point.get(year()), 1970);
            TEST_EQ(time_point.get(extended_year()), 1970);
            if(backend_name == "icu") {
                time_point = extended_year(-3);
                TEST_EQ(time_point.get(era()), 0);
                TEST_EQ(time_point.get(year()), 4);
            }

            time_point = tp_5_feb_1970_153313;
            TEST_EQ(time_point.get(month()), 1);
            TEST_EQ(time_point.get(day()), 5);
            TEST_EQ(time_point.get(day_of_year()), 36);
            TEST_EQ(time_point.get(day_of_week()), 5);
            TEST_EQ(time_point.get(day_of_week_in_month()), 1);
            time_point = date_time(a_datetime, calendar(g("ru_RU.UTF-8")));
            TEST_EQ(time_point.get(day_of_week_local()), 4);
            time_point = year(2026) + january() + day(1);
            TEST_EQ(time_point.get(day_of_week()), 5);
            TEST_EQ(time_point.get(week_of_year()), 1);
            TEST_EQ(time_point.get(week_of_month()), 1);
            time_point = day_of_week() * 1;
            TEST_EQ(time_point.get(day()), 4);
            TEST_EQ(time_point.get(week_of_year()), 1);
            TEST_EQ(time_point.get(week_of_month()), 1);
            time_point += day() * 1;
            TEST_EQ(time_point.get(week_of_year()), 2);
            TEST_EQ(time_point.get(week_of_month()), 2);

            time_point = february() + day() * 2;

            TEST_EQ(time_point.get(week_of_year()), 6);

            // cldr changes
#if BOOST_LOCALE_ICU_VERSION >= 408 && BOOST_LOCALE_ICU_VERSION <= 6000
            const bool ICU_cldr_issue = backend_name == "icu";
#else
            const bool ICU_cldr_issue = false;
#endif
            BOOST_LOCALE_START_CONST_CONDITION

            TEST_EQ(time_point.get(week_of_month()), ICU_cldr_issue ? 2 : 1);

            time_point = year(2010) + january() + day() * 3;

            TEST_EQ(time_point.get(week_of_year()), ICU_cldr_issue ? 1 : 53);

            time_point = year() * 2010 + january() + day() * 4;

            TEST_EQ(time_point.get(week_of_year()), ICU_cldr_issue ? 2 : 1);

            time_point = year() * 2010 + january() + day() * 10;

            TEST_EQ(time_point.get(week_of_year()), ICU_cldr_issue ? 2 : 1);

            time_point = year() * 2010 + january() + day() * 11;

            TEST_EQ(time_point.get(week_of_year()), ICU_cldr_issue ? 3 : 2);

            BOOST_LOCALE_END_CONST_CONDITION

            time_point = date_time(a_datetime);
            TEST_EQ(time_point.get(hour()), 15);
            TEST_EQ(date_time(a_datetime, calendar("GMT+01:00")).get(hour()), 16);
            TEST_EQ(time_point.get(hour_12()), 3);
            TEST_EQ(time_point.get(am_pm()), 1);
            TEST_EQ(time_point.get(minute()), 33);
            TEST_EQ(time_point.get(second()), 13);
            TEST_EQ(date_time(year() * 1984 + february() + day()).get(week_of_year()), 5);
            TEST_EQ(time_point.get(week_of_month()), 1);

            time_point.time(24 * 3600. * 2);

            time_point = year() * 2011;
            time_point = march();
            time_point = day() * 29;

            TEST_EQ(time_point.get(year()), 2011);
            TEST_EQ(time_point.get(month()), 2); // march
            TEST_EQ(time_point.get(day()), 29);

            date_time tp_29_march_2011 = time_point;

            time_point = year() * 2011;
            time_point = february();
            time_point = day() * 5;

            TEST_EQ(time_point.get(year()), 2011);
            TEST_EQ(time_point.get(month()), 2); // march
            TEST_EQ(time_point.get(day()), 5);

            time_point = tp_29_march_2011;

            time_point = year() * 2011 + february() + day() * 5;
            TEST_EQ(time_point.get(year()), 2011);
            TEST_EQ(time_point.get(month()), 1); // february
            TEST_EQ(time_point.get(day()), 5);

            // Difference
            TEST_EQ(time_point.difference(time_point + second(3), second()), 3);
            TEST_EQ(time_point.difference(time_point - minute(5), minute()), -5);
            TEST_EQ(time_point.difference(time_point - minute(5), second()), -5 * 60);
            TEST_EQ(time_point.difference(time_point + minute(5) - hour(3), hour()), -2);
            TEST_EQ(time_point.difference(time_point + day(42) - hour(3), day()), 41);
            TEST_EQ(time_point.difference(time_point + day(7 * 13) - hour(3), week_of_year()), 12);
            TEST_EQ(time_point.difference(time_point + day(456), day()), 456);
            TEST_EQ(time_point.difference(time_point + day(456), year()), 1);
            // Same for subtracting timepoints, i.e. syntactic sugar for the above
            TEST_EQ(((time_point + second(3)) - time_point) / second(), 3);
            TEST_EQ(((time_point - minute(5)) - time_point) / minute(), -5);
            TEST_EQ(((time_point - minute(5)) - time_point) / second(), -5 * 60);
            TEST_EQ(((time_point + minute(5) - hour(3)) - time_point) / hour(), -2);
            TEST_EQ(((time_point + day(42) - hour(3)) - time_point) / day(), 41);
            TEST_EQ(((time_point + day(7 * 13) - hour(3)) - time_point) / week_of_year(), 12);
            TEST_EQ(((time_point + day(456)) - time_point) / day(), 456);
            TEST_EQ(((time_point + day(456)) - time_point) / year(), 1);

            TEST_EQ((time_point + 2 * hour() - time_point) / minute(), 120);
            TEST_EQ((time_point + month() - time_point) / day(), 28);
            TEST_EQ((time_point + 2 * month() - (time_point + month())) / day(), 31);
            TEST_EQ((time_point + month(2) + day(3) - time_point) / month(), 2);
            TEST_EQ(day(time_point + 2 * month() - (time_point + month())), 31);
            TEST_EQ((time_point + year() * 1 - hour() * 1 - time_point) / year(), 0);
            TEST_EQ((time_point + year() * 1 - time_point) / year(), 1);
            TEST_EQ((time_point + year() * 1 + hour() * 1 - time_point) / year(), 1);
            TEST_EQ((time_point - year() * 1 + hour() * 1 - time_point) / year(), 0);
            TEST_EQ((time_point - year() * 1 - time_point) / year(), -1);
            TEST_EQ((time_point - year() * 1 - hour() * 1 - time_point) / year(), -1);
            TEST_EQ((time_point - tp_29_march_2011) / era(), 0);
            const date_time tp_morning = time_point = hour(5) + minute(7) + second(42);
            TEST_EQ(((tp_morning + am()) - tp_morning) / am_pm(), 0);
            TEST_EQ(((tp_morning + pm()) - tp_morning) / am_pm(), 1);
            // Same point
            TEST_EQ((time_point - time_point) / year(), 0);
            TEST_EQ((time_point - time_point) / month(), 0);
            TEST_EQ((time_point - time_point) / day(), 0);
            TEST_EQ((time_point - time_point) / hour(), 0);
            TEST_EQ((time_point - time_point) / minute(), 0);
            TEST_EQ((time_point - time_point) / second(), 0);
        }
        // Default constructed time_point
        {
            const time_t current_time = std::time(nullptr);
            date_time time_point_default;
            // Defaults to current time, i.e. different than a date in 1970
            date_time time_point_1970 = year(1970) + february() + day(5);
            TEST(time_point_default != time_point_1970);
            // We can not check an exact time as we can't know at which exact time the time point was recorded. So
            // only check that it refers to the same hour
            const double time_point_time = time_point_default.time();
            TEST_GE(time_point_time, current_time);
            constexpr double secsPerHour = 60 * 60;
            TEST_LE(time_point_time - current_time, secsPerHour);
            // However at least the date should match
            const tm current_time_gmt = *gmtime_wrap(&current_time);
            TEST_EQ(time_point_default.get(year()), current_time_gmt.tm_year + 1900);
            TEST_EQ(time_point_default.get(month()), current_time_gmt.tm_mon);
            TEST_EQ(time_point_default.get(day()), current_time_gmt.tm_mday);

            // Uses the current global timezone
            time_zone::global("GMT");
            date_time tp_gmt;
            time_zone::global("GMT+01:00");
            date_time tp_gmt1;
            // Both refer to the same point in time (i.e. comparison ignores timezones)
            // Unless the system clock resolution is high enough to detect that the 2 instances
            // are not created in the exact same second
            TEST((tp_gmt == tp_gmt1) || (tp_gmt1 - tp_gmt) / second() < 5);

            // But getting the hour shows the difference of 1 hour
            const int gmt_h = tp_gmt.get(hour());
            // Handle overflow to next day
            const int expected_gmt1_h = (gmt_h == tp_gmt.maximum(hour())) ? tp_gmt.minimum(hour()) : gmt_h + 1;
            TEST_EQ(expected_gmt1_h, tp_gmt1.get(hour()));
            // Adding the hour automatically handles the overflow, so this works too
            tp_gmt += hour();
            TEST_EQ(tp_gmt.get(hour()), tp_gmt1.get(hour()));
        }
    } // for loop
}

// boostinspect:noascii
