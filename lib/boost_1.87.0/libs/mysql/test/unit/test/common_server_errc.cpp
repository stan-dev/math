//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/common_server_errc.hpp>

#include <boost/test/unit_test.hpp>

#include <limits>

using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_SUITE(test_common_server_errc)

BOOST_AUTO_TEST_CASE(error_to_string_regular)
{
    struct
    {
        const char* name;
        int err;
        const char* expected_msg;
    } test_cases[] = {
        {"int_min",                  (std::numeric_limits<int>::min)(),                     "<unknown server error>"     },
        {"zero",                     0,                                                     "<unknown server error>"     },
        {"lt_min",                   999,                                                   "<unknown server error>"     },
        {"min",                      1000,                                                  "er_hashchk"                 },
        {"unused_intermediate_code", 1150,                                                  "<unknown server error>"     },
        {"db_specific",              1076,                                                  "<unknown server error>"     },
        {"regular",                  static_cast<int>(common_server_errc::er_bad_db_error), "er_bad_db_error"            },
        {"max",                      1879,                                                  "er_innodb_ft_aux_not_hex_id"},
        {"gt_max",                   1880,                                                  "<unknown server error>"     },
        {"3000_range",               3000,                                                  "<unknown server error>"     },
        {"int_max",                  (std::numeric_limits<int>::max)(),                     "<unknown server error>"     }
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            BOOST_TEST(error_code(static_cast<common_server_errc>(tc.err)).message() == tc.expected_msg);
        }
    }
}

BOOST_AUTO_TEST_CASE(error_to_string_coverage)
{
    // Check that no value causes problems.
    for (int i = 1000; i <= 1880; ++i)
    {
        BOOST_CHECK_NO_THROW(error_code(static_cast<common_server_errc>(i)).message());
    }
}

BOOST_AUTO_TEST_CASE(make_error_code)
{
    error_code code(common_server_errc::er_bad_db_error);
    BOOST_TEST(code.value() == static_cast<int>(common_server_errc::er_bad_db_error));
    BOOST_TEST(&code.category() == &get_common_server_category());  // categories are not printable
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
