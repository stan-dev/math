//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/error_categories.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/mysql_server_errc.hpp>

#include <boost/test/unit_test.hpp>

#include <limits>

using namespace boost::mysql;

namespace {

BOOST_AUTO_TEST_SUITE(test_mysql_server_errc)

std::string error_to_string(int val)
{
    error_code ec{val, get_mysql_server_category()};
    return ec.message();
}

BOOST_AUTO_TEST_CASE(error_to_string_regular)
{
    struct
    {
        const char* name;
        int code;
        const char* expected_msg;
    } test_cases[] = {
        {"int_min",           (std::numeric_limits<int>::min)(), "<unknown MySQL-specific server error>"  },
        {"zero",              0,                                 "<unknown MySQL-specific server error>"  },
        {"common_min",        1000,                              "<unknown MySQL-specific server error>"  },
        {"common_repurposed", 1076,                              "er_ready"                               }, // != in MariaDB
        {"common_max",        1879,                              "<unknown MySQL-specific server error>"  },
        {"specific1_min",     1880,                              "er_old_temporals_upgraded"              },
        {"specific1_max",     1888,                              "er_found_missing_gtids"                 },
        {"specific1_gtmax",   1889,                              "<unknown MySQL-specific server error>"  },
        {"specific2_ltmin",   2999,                              "<unknown MySQL-specific server error>"  },
        {"specific2_min",     3000,                              "er_file_corrupt"                        },
        {"specific2_regular", 3015,                              "er_engine_out_of_memory"                },
        {"specific2_max",     3238,                              "er_aes_invalid_kdf_option_size"         },
        {"specific2_gtmax",   3239,                              "<unknown MySQL-specific server error>"  },
        {"specific3_ltmin",   3499,                              "<unknown MySQL-specific server error>"  },
        {"specific3_min",     3500,                              "er_unsupport_compressed_temporary_table"},
        {"specific3_regular", 3532,                              "er_rename_role"                         },
        {"gt_max",            5000,                              "<unknown MySQL-specific server error>"  },
        {"int_max",           (std::numeric_limits<int>::max)(), "<unknown MySQL-specific server error>"  },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name) { BOOST_TEST(error_to_string(tc.code) == tc.expected_msg); }
    }
}

BOOST_AUTO_TEST_CASE(error_to_string_coverage)
{
    // Check that no value causes problems.
    // Ensure that all branches of the switch/case are covered
    // Valid error ranges are 1000-2000 and 3000-5000
    for (int i = 1000; i < 2000; ++i)
    {
        BOOST_CHECK_NO_THROW(error_to_string(i));
    }
    for (int i = 3000; i < 5000; ++i)
    {
        BOOST_CHECK_NO_THROW(error_to_string(i));
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace
