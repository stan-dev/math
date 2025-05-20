//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/blob.hpp>
#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/date.hpp>
#include <boost/mysql/datetime.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/time.hpp>

#include <boost/config.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/optional/optional.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <limits>
#include <string>

#include "format_common.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"

#ifdef __cpp_lib_string_view
#include <string_view>
#endif
#ifdef __cpp_lib_optional
#include <optional>
#endif

using namespace boost::mysql;
using namespace boost::mysql::test;
using mysql_time = boost::mysql::time;

//
// Verify that formatting individual values work. This is tested using format_sql
// because it's convenient, but also covers basic_format_context
//
BOOST_AUTO_TEST_SUITE(test_format_sql_individual_value)

constexpr format_options opts{utf8mb4_charset, true};
constexpr const char* single_fmt = "SELECT {};";
constexpr const char* identifier_fmt = "SELECT {:i} FROM myt";
constexpr const char* raw_fmt = "SELECT {:r};";

BOOST_AUTO_TEST_CASE(null_)
{
    // nullptr interpreted as NULL
    BOOST_TEST(format_sql(opts, single_fmt, nullptr) == "SELECT NULL;");
}

BOOST_AUTO_TEST_CASE(signed_char)
{
    BOOST_TEST(format_sql(opts, single_fmt, (signed char)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (signed char)-1) == "SELECT -1;");
    BOOST_TEST(format_sql(opts, single_fmt, (signed char)-128) == "SELECT -128;");
    BOOST_TEST(format_sql(opts, single_fmt, (signed char)127) == "SELECT 127;");
}

BOOST_AUTO_TEST_CASE(unsigned_char)
{
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned char)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned char)0) == "SELECT 0;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned char)0xff) == "SELECT 255;");
}

BOOST_AUTO_TEST_CASE(short_)
{
    BOOST_TEST(format_sql(opts, single_fmt, (short)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (short)-1) == "SELECT -1;");
    BOOST_TEST(format_sql(opts, single_fmt, (short)-32768) == "SELECT -32768;");
    BOOST_TEST(format_sql(opts, single_fmt, (short)32767) == "SELECT 32767;");
}

BOOST_AUTO_TEST_CASE(unsigned_short)
{
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned short)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned short)0) == "SELECT 0;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned short)0xffff) == "SELECT 65535;");
}

BOOST_AUTO_TEST_CASE(int_)
{
    BOOST_TEST(format_sql(opts, single_fmt, (int)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (int)-1) == "SELECT -1;");
    BOOST_TEST(format_sql(opts, single_fmt, (int)-0x7fffffff - 1) == "SELECT -2147483648;");
    BOOST_TEST(format_sql(opts, single_fmt, (int)0x7fffffff) == "SELECT 2147483647;");
}

BOOST_AUTO_TEST_CASE(unsigned_int)
{
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned int)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned int)0) == "SELECT 0;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned int)0xffffffff) == "SELECT 4294967295;");
}

BOOST_AUTO_TEST_CASE(long_)
{
    BOOST_TEST(format_sql(opts, single_fmt, (long)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (long)-1) == "SELECT -1;");
    BOOST_TEST(format_sql(opts, single_fmt, (long)0x7fffffff) == "SELECT 2147483647;");
}

BOOST_AUTO_TEST_CASE(unsigned_long)
{
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned long)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned long)0) == "SELECT 0;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned long)0xffffffff) == "SELECT 4294967295;");
}

BOOST_AUTO_TEST_CASE(long_long)
{
    BOOST_TEST(format_sql(opts, single_fmt, (long long)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (long long)-1) == "SELECT -1;");
    BOOST_TEST(
        format_sql(opts, single_fmt, (long long)(-0x7fffffffffffffff - 1)) == "SELECT -9223372036854775808;"
    );
    BOOST_TEST(format_sql(opts, single_fmt, (long long)0x7fffffffffffffff) == "SELECT 9223372036854775807;");
}

BOOST_AUTO_TEST_CASE(unsigned_long_long)
{
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned long long)42) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, (unsigned long long)0) == "SELECT 0;");
    BOOST_TEST(
        format_sql(opts, single_fmt, (unsigned long long)0xffffffffffffffff) == "SELECT 18446744073709551615;"
    );
}

BOOST_AUTO_TEST_CASE(bool_)
{
    BOOST_TEST(format_sql(opts, single_fmt, true) == "SELECT 1;");
    BOOST_TEST(format_sql(opts, single_fmt, false) == "SELECT 0;");
}

BOOST_AUTO_TEST_CASE(double_)
{
    // doubles have many different cases that may cause trouble
    struct
    {
        string_view name;
        double value;
        string_view expected;
    } test_cases[] = {
        {"regular",               4.2,                      "4.2e+00"                 },
        {"regular_precision",     4.298238239237823287327,  "4.298238239237823e+00"   },
        {"exp",                   5.1e+23,                  "5.1e+23"                 },
        {"exp_precision",         4.2982382392378232e+67,   "4.2982382392378234e+67"  },
        {"max",                   1.7976931348623157e+308,  "1.7976931348623157e+308" },
        {"regular_neg",           -4.2,                     "-4.2e+00"                },
        {"regular_precision_neg", -4.298238239237823287327, "-4.298238239237823e+00"  },
        {"exp_neg",               -5.1e+23,                 "-5.1e+23"                },
        {"max_neg",               -1.7976931348623157e+308, "-1.7976931348623157e+308"},
        {"zero",                  0.0,                      "0e+00"                   },
        {"zero_neg",              -0.0,                     "-0e+00"                  },
        {"expneg",                4.2e-12,                  "4.2e-12"                 },
        {"expneg_precision",      4.2872383293922839e-45,   "4.2872383293922836e-45"  },
        {"min",                   2.2250738585072014e-308,  "2.2250738585072014e-308" },
        {"min_neg",               -2.2250738585072014e-308, "-2.2250738585072014e-308"},
        {"denorm",                -4.2872383293922839e-309, "-4.287238329392283e-309" },
        {"min_denorm",            5e-324,                   "5e-324"                  },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            auto str = format_sql(opts, "{}", tc.value);
            BOOST_TEST(str == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(float_)
{
    // floats are converted to double before formatting, since MySQL
    // interprets all floating point literals as doubles
    struct
    {
        string_view name;
        float value;
        string_view expected;
    } test_cases[] = {
        {"regular",               4.2f,                      "4.199999809265137e+00"  },
        {"regular_precision",     4.298238239237823287327f,  "4.298238277435303e+00"  },
        {"exp",                   5.1e+23f,                  "5.100000157096095e+23"  },
        {"exp_precision",         4.2982382392378232e+35f,   "4.298238339685548e+35"  },
        {"max",                   3.4028234663852886e+38f,   "3.4028234663852886e+38" },
        {"regular_neg",           -4.2f,                     "-4.199999809265137e+00" },
        {"regular_precision_neg", -4.298238239237823287327f, "-4.298238277435303e+00" },
        {"exp_neg",               -5.1e+23f,                 "-5.100000157096095e+23" },
        {"max_neg",               -3.4028234663852886e+38f,  "-3.4028234663852886e+38"},
        {"zero",                  0.0f,                      "0e+00"                  },
        {"zero_neg",              -0.0f,                     "-0e+00"                 },
        {"expneg",                4.2e-12f,                  "4.200000156689976e-12"  },
        {"expneg_precision",      4.2872383293922839e-23f,   "4.2872384543670994e-23" },
        {"min",                   1.1754944e-38f,            "1.1754943508222875e-38" },
        {"min_neg",               -1.1754944e-38f,           "-1.1754943508222875e-38"},
        {"denorm",                -4.2872383293922839e-39f,  "-4.287239020438634e-39" },
        {"min_denorm",            1.401298464324817e-45f,    "1.401298464324817e-45"  },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            auto str = format_sql(opts, "{}", tc.value);
            BOOST_TEST(str == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(string_literal)
{
    // As values
    BOOST_TEST(format_sql(opts, single_fmt, "abc") == "SELECT 'abc';");
    BOOST_TEST(format_sql(opts, single_fmt, "abc'\\ OR 1=1") == "SELECT 'abc\\'\\\\ OR 1=1';");
    BOOST_TEST(format_sql(opts, single_fmt, "hola \xc3\xb1!") == "SELECT 'hola \xc3\xb1!';");
    BOOST_TEST(format_sql(opts, single_fmt, "") == "SELECT '';");

    // As identifiers
    BOOST_TEST(format_sql(opts, identifier_fmt, "myfield") == "SELECT `myfield` FROM myt");
    BOOST_TEST(format_sql(opts, identifier_fmt, "inj`ect'ion") == "SELECT `inj``ect'ion` FROM myt");
    BOOST_TEST(
        format_sql(opts, identifier_fmt, "mo`e\\inj``ectionatt\nemmpts`") ==
        "SELECT `mo``e\\inj````ectionatt\nemmpts``` FROM myt"
    );

    // Empty identifiers are not valid in MySQL but they shouldn't case formatting problems.
    // They are correctly rejected by MySQL (they don't cause problems)
    BOOST_TEST(format_sql(opts, identifier_fmt, "") == "SELECT `` FROM myt");

    // As raw
    BOOST_TEST(format_sql(opts, raw_fmt, "abc") == "SELECT abc;");              // regular
    BOOST_TEST(format_sql(opts, raw_fmt, "") == "SELECT ;");                    // empty
    BOOST_TEST(format_sql(opts, raw_fmt, "a\\'\"b`c") == "SELECT a\\'\"b`c;");  // we don't escape
    BOOST_TEST(format_sql(opts, raw_fmt, "a\xff bc") == "SELECT a\xff bc;");    // we don't check charset
}

BOOST_AUTO_TEST_CASE(c_str)
{
    BOOST_TEST(format_sql(opts, single_fmt, static_cast<const char*>("abc")) == "SELECT 'abc';");
    BOOST_TEST(format_sql(opts, single_fmt, static_cast<const char*>("")) == "SELECT '';");
}

BOOST_AUTO_TEST_CASE(string)
{
    std::string lval = "I'm an lvalue";
    const std::string clval = "I'm const";
    BOOST_TEST(format_sql(opts, single_fmt, lval) == "SELECT 'I\\'m an lvalue';");
    BOOST_TEST(format_sql(opts, single_fmt, clval) == "SELECT 'I\\'m const';");
    BOOST_TEST(format_sql(opts, single_fmt, std::string("abc")) == "SELECT 'abc';");
    BOOST_TEST(format_sql(opts, single_fmt, std::string()) == "SELECT '';");
    BOOST_TEST(format_sql(opts, single_fmt, string_with_alloc("abc'")) == "SELECT 'abc\\'';");

    // specifiers work
    BOOST_TEST(format_sql(opts, identifier_fmt, lval) == "SELECT `I'm an lvalue` FROM myt");
    BOOST_TEST(format_sql(opts, identifier_fmt, clval) == "SELECT `I'm const` FROM myt");
    BOOST_TEST(format_sql(opts, identifier_fmt, std::string("abc")) == "SELECT `abc` FROM myt");
    BOOST_TEST(format_sql(opts, identifier_fmt, string_with_alloc("abc")) == "SELECT `abc` FROM myt");
}

BOOST_AUTO_TEST_CASE(string_view_)
{
    BOOST_TEST(format_sql(opts, single_fmt, string_view("abc")) == "SELECT 'abc';");
    BOOST_TEST(format_sql(opts, single_fmt, string_view("abc'\\ OR 1=1")) == "SELECT 'abc\\'\\\\ OR 1=1';");
    BOOST_TEST(format_sql(opts, single_fmt, string_view()) == "SELECT '';");

    // specifiers work
    BOOST_TEST(format_sql(opts, identifier_fmt, string_view("abc")) == "SELECT `abc` FROM myt");
}

#ifdef __cpp_lib_string_view
BOOST_AUTO_TEST_CASE(std_string_view)
{
    BOOST_TEST(format_sql(opts, single_fmt, std::string_view("abc")) == "SELECT 'abc';");
    BOOST_TEST(
        format_sql(opts, single_fmt, std::string_view("abc'\\ OR 1=1")) == "SELECT 'abc\\'\\\\ OR 1=1';"
    );
    BOOST_TEST(format_sql(opts, single_fmt, std::string_view()) == "SELECT '';");

    // specifiers work
    BOOST_TEST(format_sql(opts, identifier_fmt, std::string_view("abc")) == "SELECT `abc` FROM myt");
}
#endif

// blob: encoded as a hex string
BOOST_AUTO_TEST_CASE(blob_)
{
    blob lval{0x01, 0x00, 0x5c};
    const blob clval{0x20, 0x71, 0xff};
    BOOST_TEST(format_sql(opts, single_fmt, lval) == "SELECT x'01005c';");
    BOOST_TEST(format_sql(opts, single_fmt, clval) == "SELECT x'2071ff';");
    BOOST_TEST(format_sql(opts, single_fmt, blob{0x00, 0x2c}) == "SELECT x'002c';");
}

BOOST_AUTO_TEST_CASE(blob_coverage)
{
    struct
    {
        string_view name;
        const blob input;
        string_view expected;
    } test_cases[] = {
        {"empty",           {},                                                                            "SELECT x'';"        },
        {"injection_chars", {0x5c, 0x5c, 0x27},                                                            "SELECT x'5c5c27';"  }, // 5c = backslash, 27 = single quote
        {"all_zeros",       {0x00, 0x00, 0x00, 0x00},                                                      "SELECT x'00000000';"},

        // Check that we encode all possible byte values correctly
        {"bytes_00_3f",
         {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
          0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
          0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
          0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f},
         "SELECT x'000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
         "202122232425262728292a2b2c2d2e2f303132333435363738393a3b3c3d3e3f';"                                                   },
        {"bytes_40_7f",
         {0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f,
          0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f,
          0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f,
          0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f},
         "SELECT x'404142434445464748494a4b4c4d4e4f505152535455565758595a5b5c5d5e5f606162"
         "636465666768696a6b6c6d6e6f707172737475767778797a7b7c7d7e7f';"                                                         },
        {"bytes_80_bf",
         {0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
          0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
          0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
          0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf},
         "SELECT x'808182838485868788898a8b8c8d8e8f909192939495969798999a9b9c9d9e9f"
         "a0a1a2a3a4a5a6a7a8a9aaabacadaeafb0b1b2b3b4b5b6b7b8b9babbbcbdbebf';"                                                   },
        {"bytes_c0_ff",
         {0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
          0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
          0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
          0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff},
         "SELECT x'c0c1c2c3c4c5c6c7c8c9cacbcccdcecfd0d1d2d3d4d5d6d7d8d9dadbdcdddedf"
         "e0e1e2e3e4e5e6e7e8e9eaebecedeeeff0f1f2f3f4f5f6f7f8f9fafbfcfdfeff';"                                                   },

        // We use a 64 byte buffer for the formatting operation. Update these if the buffer size changes.
        {"31_bytes",
         {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
          0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f},
         "SELECT x'000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f';"                                          },
        {"32_bytes",
         {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
          0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
          0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20},
         "SELECT x'000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20';"                                        },
        {"33_bytes",
         {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b,
          0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
          0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21},
         "SELECT x'000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f2021';"                                      },
        {"63_bytes",
         {0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
          0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
          0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
          0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe},
         "SELECT x'808182838485868788898a8b8c8d8e8f909192939495969798999a9b9c9d9e9f"
         "a0a1a2a3a4a5a6a7a8a9aaabacadaeafb0b1b2b3b4b5b6b7b8b9babbbcbdbe';"                                                     },
        {"64_bytes",
         {0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
          0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
          0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
          0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf},
         "SELECT x'808182838485868788898a8b8c8d8e8f909192939495969798999a9b9c9d9e9f"
         "a0a1a2a3a4a5a6a7a8a9aaabacadaeafb0b1b2b3b4b5b6b7b8b9babbbcbdbebf';"                                                   },
        {"65_bytes",
         {0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c,
          0x8d, 0x8e, 0x8f, 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99,
          0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f, 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6,
          0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf, 0xb0, 0xb1, 0xb2, 0xb3,
          0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf, 0xc0},
         "SELECT x'808182838485868788898a8b8c8d8e8f909192939495969798999a9b9c9d9e9f"
         "a0a1a2a3a4a5a6a7a8a9aaabacadaeafb0b1b2b3b4b5b6b7b8b9babbbcbdbebfc0';"                                                 },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name) { BOOST_TEST(format_sql(opts, single_fmt, tc.input) == tc.expected); }
    }
}

BOOST_AUTO_TEST_CASE(blob_view_)
{
    BOOST_TEST(format_sql(opts, single_fmt, makebv("hello\\")) == "SELECT x'68656c6c6f5c';");
    BOOST_TEST(format_sql(opts, single_fmt, makebv("hello \xc3\xb1!")) == "SELECT x'68656c6c6f20c3b121';");
    BOOST_TEST(format_sql(opts, single_fmt, makebv("hello \xc3'!")) == "SELECT x'68656c6c6f20c32721';");
    BOOST_TEST(format_sql(opts, single_fmt, blob_view()) == "SELECT x'';");
}

BOOST_AUTO_TEST_CASE(blob_array)
{
    // Collections of unsigned char are formatted as blob if they're convertible to span<const unsigned char>
    std::array<unsigned char, 4> arr{
        {5, 1, 0, 0xab}
    };
    const std::array<unsigned char, 4> carr{
        {0xde, 0xad, 0xbe, 0xef}
    };

    BOOST_TEST(format_sql(opts, single_fmt, arr) == "SELECT x'050100ab';");
    BOOST_TEST(format_sql(opts, single_fmt, carr) == "SELECT x'deadbeef';");
}

BOOST_AUTO_TEST_CASE(date_)
{
    BOOST_TEST(format_sql(opts, single_fmt, date(2021, 1, 20)) == "SELECT '2021-01-20';");
    BOOST_TEST(format_sql(opts, single_fmt, date()) == "SELECT '0000-00-00';");
    BOOST_TEST(format_sql(opts, single_fmt, date(0xffff, 0xff, 0xff)) == "SELECT '65535-255-255';");
}

BOOST_AUTO_TEST_CASE(datetime_)
{
    BOOST_TEST(format_sql(opts, single_fmt, datetime(2021, 1, 20)) == "SELECT '2021-01-20 00:00:00.000000';");
    BOOST_TEST(
        format_sql(opts, single_fmt, datetime(1998, 1, 1, 21, 3, 5, 12)) ==
        "SELECT '1998-01-01 21:03:05.000012';"
    );
    BOOST_TEST(format_sql(opts, single_fmt, datetime()) == "SELECT '0000-00-00 00:00:00.000000';");
    BOOST_TEST(
        format_sql(opts, single_fmt, datetime(0xffff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xffffffff)) ==
        "SELECT '65535-255-255 255:255:255.4294967295';"
    );
}

BOOST_AUTO_TEST_CASE(time_)
{
    BOOST_TEST(format_sql(opts, single_fmt, maket(127, 1, 10, 123)) == "SELECT '127:01:10.000123';");
    BOOST_TEST(format_sql(opts, single_fmt, -maket(9, 1, 10)) == "SELECT '-09:01:10.000000';");
    BOOST_TEST(format_sql(opts, single_fmt, mysql_time()) == "SELECT '00:00:00.000000';");
    BOOST_TEST(format_sql(opts, single_fmt, (mysql_time::min)()) == "SELECT '-2562047788:00:54.775808';");
    BOOST_TEST(format_sql(opts, single_fmt, (mysql_time::max)()) == "SELECT '2562047788:00:54.775807';");
}

BOOST_AUTO_TEST_CASE(duration)
{
    // Durations work as long as they're compatible with time
    using namespace std::chrono;

    BOOST_TEST(format_sql(opts, single_fmt, hours(21)) == "SELECT '21:00:00.000000';");
    BOOST_TEST(format_sql(opts, single_fmt, minutes(3)) == "SELECT '00:03:00.000000';");
    BOOST_TEST(format_sql(opts, single_fmt, seconds(-10)) == "SELECT '-00:00:10.000000';");
    BOOST_TEST(format_sql(opts, single_fmt, milliseconds(-9)) == "SELECT '-00:00:00.009000';");
    BOOST_TEST(format_sql(opts, single_fmt, microseconds(3214)) == "SELECT '00:00:00.003214';");
}

BOOST_AUTO_TEST_CASE(field_view_)
{
    field referenced("def\\");
    BOOST_TEST(format_sql(opts, single_fmt, field_view()) == "SELECT NULL;");
    BOOST_TEST(format_sql(opts, single_fmt, field_view(42)) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, field_view("'abc'")) == "SELECT '\\'abc\\'';");
    BOOST_TEST(format_sql(opts, single_fmt, field_view(referenced)) == "SELECT 'def\\\\';");
}

BOOST_AUTO_TEST_CASE(field_)
{
    field f_lval("hol\"a");
    const field f_clval(42);
    BOOST_TEST(format_sql(opts, single_fmt, field()) == "SELECT NULL;");
    BOOST_TEST(format_sql(opts, single_fmt, field(4.2)) == "SELECT 4.2e+00;");
    BOOST_TEST(format_sql(opts, single_fmt, f_lval) == "SELECT 'hol\\\"a';");
    BOOST_TEST(format_sql(opts, single_fmt, f_clval) == "SELECT 42;");
}

BOOST_AUTO_TEST_CASE(boost_optional)
{
    boost::optional<std::string> o_lval("abc");
    boost::optional<const std::string> co_lval("ab'c");
    const boost::optional<std::string> o_clval("\\");
    const boost::optional<const std::string> co_clval("abdef");
    BOOST_TEST(format_sql(opts, single_fmt, boost::optional<int>()) == "SELECT NULL;");
    BOOST_TEST(format_sql(opts, single_fmt, boost::optional<int>(42)) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, o_lval) == "SELECT 'abc';");
    BOOST_TEST(format_sql(opts, single_fmt, co_lval) == "SELECT 'ab\\'c';");
    BOOST_TEST(format_sql(opts, single_fmt, o_clval) == "SELECT '\\\\';");
    BOOST_TEST(format_sql(opts, single_fmt, co_clval) == "SELECT 'abdef';");
}

#ifdef __cpp_lib_optional
BOOST_AUTO_TEST_CASE(std_optional)
{
    std::optional<std::string> o_lval("abc");
    std::optional<const std::string> co_lval("ab'c");
    const std::optional<std::string> o_clval("\\");
    const std::optional<const std::string> co_clval("abdef");
    BOOST_TEST(format_sql(opts, single_fmt, std::optional<int>()) == "SELECT NULL;");
    BOOST_TEST(format_sql(opts, single_fmt, std::optional<int>(42)) == "SELECT 42;");
    BOOST_TEST(format_sql(opts, single_fmt, o_lval) == "SELECT 'abc';");
    BOOST_TEST(format_sql(opts, single_fmt, co_lval) == "SELECT 'ab\\'c';");
    BOOST_TEST(format_sql(opts, single_fmt, o_clval) == "SELECT '\\\\';");
    BOOST_TEST(format_sql(opts, single_fmt, co_clval) == "SELECT 'abdef';");
}
#endif

//
// Errors when formatting individual fields
//

// Make code less verbose
constexpr auto spec_err = client_errc::format_string_invalid_specifier;

BOOST_AUTO_TEST_CASE(null_error)
{
    // Specifiers rejected
    BOOST_TEST(format_single_error("SELECT {:i}", nullptr) == spec_err);
}

BOOST_AUTO_TEST_CASE(integers_error)
{
    // Specifiers rejected
    BOOST_TEST(format_single_error("SELECT {:i}", (signed char)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (unsigned char)0xff) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (short)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:r}", (unsigned short)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (int)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:d}", (int)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {: }", (int)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {::}", (int)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (unsigned int)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (long)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (unsigned long)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (long long)42) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i}", (unsigned long long)42) == spec_err);
}

BOOST_AUTO_TEST_CASE(bool_error)
{
    // Specifiers rejected
    BOOST_TEST(format_single_error("SELECT {:i}", true) == spec_err);
    BOOST_TEST(format_single_error("SELECT {:r}", false) == spec_err);
}

BOOST_AUTO_TEST_CASE(double_error)
{
    using double_limits = std::numeric_limits<double>;

    // double inf and nan not supported
    BOOST_TEST(format_single_error("{}", double_limits::infinity()) == client_errc::unformattable_value);
    BOOST_TEST(format_single_error("{}", -double_limits::infinity()) == client_errc::unformattable_value);
    BOOST_TEST(format_single_error("{}", double_limits::quiet_NaN()) == client_errc::unformattable_value);

    // Specifiers rejected
    BOOST_TEST(format_single_error("SELECT {:i}", 21.0e10) == spec_err);
    BOOST_TEST(format_single_error("SELECT {: f}", 0.0) == spec_err);
}

BOOST_AUTO_TEST_CASE(float_error)
{
    using float_limits = std::numeric_limits<float>;

    // float inf and nan not supported
    BOOST_TEST(format_single_error("{}", float_limits::infinity()) == client_errc::unformattable_value);
    BOOST_TEST(format_single_error("{}", -float_limits::infinity()) == client_errc::unformattable_value);
    BOOST_TEST(format_single_error("{}", float_limits::quiet_NaN()) == client_errc::unformattable_value);

    // Specifiers rejected
    BOOST_TEST(format_single_error("SELECT {:i}", 4.2f) == spec_err);
    BOOST_TEST(format_single_error("SELECT {: f}", 4.2f) == spec_err);
}

BOOST_AUTO_TEST_CASE(string_error)
{
    // strings with invalid characters
    BOOST_TEST(format_single_error("{}", "a\xc3'") == client_errc::invalid_encoding);
    BOOST_TEST(format_single_error("{}", "a\xc3\''") == client_errc::invalid_encoding);
    BOOST_TEST(format_single_error("{}", "a\xff\xff") == client_errc::invalid_encoding);

    // identifiers with invalid characters
    BOOST_TEST(format_single_error("{:i}", "a\xd8") == client_errc::invalid_encoding);
    BOOST_TEST(format_single_error("{:i}", "a\xc3 ") == client_errc::invalid_encoding);

    // unknown specifiers are rejected
    BOOST_TEST(format_single_error("SELECT {:x}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:s}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:d}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:>}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {::}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:id}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:ir}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:ri}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:sd}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i:}", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:i }", "abc") == spec_err);
    BOOST_TEST(format_single_error("SELECT {:ivery long [value] with\" quotes'}", "abc") == spec_err);
}

BOOST_AUTO_TEST_CASE(blob_error)
{
    // blobs reject specifiers
    blob lval{0x01, 0x00, 0x5c};
    const blob clval{0x20, 0x71, 0xff};
    BOOST_TEST(format_single_error("{:i}", lval) == spec_err);
    BOOST_TEST(format_single_error("{:i}", clval) == spec_err);
    BOOST_TEST(format_single_error("{:i}", blob_view(clval)) == spec_err);
}

BOOST_AUTO_TEST_CASE(date_error)
{
    // date reject specifiers
    BOOST_TEST(format_single_error("{:i}", date(2021, 1, 20)) == spec_err);
}

BOOST_AUTO_TEST_CASE(datetime_error)
{
    // datetime reject specifiers
    BOOST_TEST(format_single_error("{:i}", datetime(1998, 1, 1, 21, 3, 5, 12)) == spec_err);
}

BOOST_AUTO_TEST_CASE(duration_error)
{
    // durations reject specifiers
    BOOST_TEST(format_single_error("{:i}", maket(9, 1, 10)) == spec_err);
    BOOST_TEST(format_single_error("{:i}", std::chrono::hours(9)) == spec_err);
}

BOOST_AUTO_TEST_CASE(field_view_error)
{
    // field_view rejects specifiers, even if the underlying type would support them
    BOOST_TEST(format_single_error("{:i}", field_view()) == spec_err);
    BOOST_TEST(format_single_error("{:i}", field_view("abc")) == spec_err);
    BOOST_TEST(format_single_error("{:r}", field_view("abc")) == spec_err);
    BOOST_TEST(format_single_error("{:i}", field_view(42)) == spec_err);

    // Errors applicable to the contained type
    BOOST_TEST(format_single_error("{}", field_view("a\xc3'")) == client_errc::invalid_encoding);
    BOOST_TEST(
        format_single_error("{}", field_view(std::numeric_limits<double>::infinity())) ==
        client_errc::unformattable_value
    );
}

BOOST_AUTO_TEST_CASE(field_error)
{
    // Same as field_view
    BOOST_TEST(format_single_error("{:i}", field("abc")) == spec_err);
    BOOST_TEST(format_single_error("{}", field("a\xc3'")) == client_errc::invalid_encoding);
    BOOST_TEST(
        format_single_error("{}", field(std::numeric_limits<double>::infinity())) ==
        client_errc::unformattable_value
    );
}

// Apply the same test to boost and std optional types
template <template <class> class Optional>
void optional_error_test()
{
    Optional<std::string> o_lval("abc");
    Optional<const std::string> co_lval("ab'c");
    const Optional<std::string> o_clval("\\");
    const Optional<const std::string> co_clval("abdef");

    // optionals rejects specifiers, even if the underlying type would support them
    BOOST_TEST(format_single_error("{:i}", o_lval) == spec_err);
    BOOST_TEST(format_single_error("{:i}", co_lval) == spec_err);
    BOOST_TEST(format_single_error("{:i}", o_clval) == spec_err);
    BOOST_TEST(format_single_error("{:i}", co_clval) == spec_err);
    BOOST_TEST(format_single_error("{:i}", Optional<std::string>("abc")) == spec_err);
    BOOST_TEST(format_single_error("{:i}", Optional<std::string>()) == spec_err);
    BOOST_TEST(format_single_error("{:r}", Optional<std::string>("abc")) == spec_err);
    BOOST_TEST(format_single_error("{:i}", Optional<int>(42)) == spec_err);

    // Errors applicable to the contained type
    BOOST_TEST(
        format_single_error("{}", Optional<std::string>("b\xff\xff")) == client_errc::invalid_encoding
    );
    BOOST_TEST(
        format_single_error("{}", Optional<double>(std::numeric_limits<double>::infinity())) ==
        client_errc::unformattable_value
    );
}

BOOST_AUTO_TEST_CASE(boost_optional_error) { optional_error_test<boost::optional>(); }

#ifdef __cpp_lib_optional
BOOST_AUTO_TEST_CASE(std_optional_error) { optional_error_test<std::optional>(); }
#endif

BOOST_AUTO_TEST_SUITE_END()
