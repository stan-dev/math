//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/escape_string.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <string>

#include "test_common/create_basic.hpp"
#include "test_common/printing.hpp"
#include "test_unit/custom_allocator.hpp"
#include "test_unit/ff_charset.hpp"

using namespace boost::mysql;

BOOST_AUTO_TEST_SUITE(test_escape_string)

//
// Escaping using backslashes
//
BOOST_AUTO_TEST_CASE(backslashes_utf8mb4_valid)
{
    char all_ascii_noescape_storage[] = {
        '\x01', '\x02', '\x03', '\x04', '\x05', '\x06', '\x07', '\x08', '\x09', '\x0b', '\x0c', '\x0e',
        '\x0f', '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17', '\x18', '\x19', '\x1b',
        '\x1c', '\x1d', '\x1e', '\x1f', ' ',    '!',    '#',    '%',    '&',    '(',    ')',    '*',
        '+',    ',',    '-',    '.',    '/',    '0',    '1',    '2',    '3',    '4',    '5',    '6',
        '7',    '8',    '9',    ':',    ';',    '<',    '=',    '>',    '?',    '@',    'A',    'B',
        'C',    'D',    'E',    'F',    'G',    'H',    'I',    'J',    'K',    'L',    'M',    'N',
        'O',    'P',    'Q',    'R',    'S',    'T',    'U',    'V',    'W',    'X',    'Y',    'Z',
        '[',    ']',    '^',    '_',    '`',    'a',    'b',    'c',    'd',    'e',    'f',    'g',
        'h',    'i',    'j',    'k',    'l',    'm',    'n',    'o',    'p',    'q',    'r',    's',
        't',    'u',    'v',    'w',    'x',    'y',    'z',    '{',    '|',    '}',    '~',
    };
    string_view all_ascii_noescape(all_ascii_noescape_storage, sizeof(all_ascii_noescape_storage));

    struct
    {
        string_view name;
        string_view input;
        string_view expected;
    } test_cases[] = {
        {"empty",                  "",                                      ""                                             },
        {"no_escape_ascii",        "this is A test string",                 "this is A test string"                        },
        {"all_noescape_ascii",     all_ascii_noescape,                      all_ascii_noescape                             },
        {"escape_dquote",          R"(this has "dquotes")",                 R"(this has \"dquotes\")"                      },
        {"escape_squote",          R"(this has 'squotes')",                 R"(this has \'squotes\')"                      },
        {"escape_backslash",       R"(this has \a backslash\)",             R"(this has \\a backslash\\)"                  },
        {"escape_null",            test::makesv("this has \0 null \0"),     R"(this has \0 null \0)"                       },
        {"escape_ctrlz",           "this has \x1a Ctrl+Z \x1a",             R"(this has \Z Ctrl+Z \Z)"                     },
        {"escape_newline",         "this has \n newline \n",                R"(this has \n newline \n)"                    },
        {"escape_carriage_return", "this has \r car ret \r",                R"(this has \r car ret \r)"                    },
        {"all_escape_chars",       "\"''\\\"\n\r\\",                        R"(\"\'\'\\\"\n\r\\)"                          },
        {"start_escape_char",      "\"abc",                                 R"(\"abc)"                                     },
        {"end_escape_char",        "abc\"",                                 R"(abc\")"                                     },
        {"single_escape_char",     "'",                                     "\\'"                                          },
        {"utf8_2byte",             "2byte \" \xc3\xb1 UTF-8\\ \xc3\xb2 \\", "2byte \\\" \xc3\xb1 UTF-8\\\\ \xc3\xb2 \\\\"
        },
        {"utf8_3byte",             "3byte '\xef\xbf\xbf UTF-8'",            "3byte \\'\xef\xbf\xbf UTF-8\\'"               },
        {"utf8_4byte",             "4byte \r'\xf0\x90\x80\x80 UTF-8\n",     "4byte \\r\\'\xf0\x90\x80\x80 UTF-8\\n"        },
 // Some typical injection payloads
        {"injection_1",            R"(\\)",                                 R"(\\\\)"                                      },
        {"injection_2",            R"(' or ")",                             R"(\' or \")"                                  },
        {"injection_3",            R"(-- or #)",                            R"(-- or #)"                                   },
        {"injection_4",            R"(' OR '1)",                            R"(\' OR \'1)"                                 },
        {"injection_5",            R"(' OR 1 -- -)",                        R"(\' OR 1 -- -)"                              },
        {"injection_6",            R"(" OR "" = ")",                        R"(\" OR \"\" = \")"                           },
        {"injection_7",            R"(" OR 1 = 1 -- -)",                    R"(\" OR 1 = 1 -- -)"                          },
        {"injection_8",            R"(' OR '' = ')",                        R"(\' OR \'\' = \')"                           },
        {"injection_9",            R"('=')",                                R"(\'=\')"                                     },
        {"injection_10",           R"('LIKE')",                             R"(\'LIKE\')"                                  },
        {"injection_11",           R"('=0--+)",                             R"(\'=0--+)"                                   },
        {"injection_12",           R"(' OR 'x'='x)",                        R"(\' OR \'x\'=\'x)"                           },
        {"injection_13",           R"(' AND id IS NULL; --)",               R"(\' AND id IS NULL; --)"                     },
        {"injection_14",           R"('''''''''''''UNION SELECT '2)",       R"(\'\'\'\'\'\'\'\'\'\'\'\'\'UNION SELECT \'2)"},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            std::string output = "abc";
            auto ec = escape_string(tc.input, {utf8mb4_charset, true}, quoting_context::double_quote, output);
            BOOST_TEST(ec == error_code());
            BOOST_TEST(output == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(backslashes_utf8mb4_invalid)
{
    // \xc0\x80 is an overlong 0 character
    std::string output = "abc";
    auto ec = escape_string(
        "This 'has' invalid \xc0\x80 chars\\",
        {utf8mb4_charset, true},
        quoting_context::double_quote,
        output
    );
    BOOST_TEST(ec == client_errc::invalid_encoding);
}

BOOST_AUTO_TEST_CASE(backslashes_multibyte_ascii_compatible_chars)
{
    // Edge cases for encodings allowing character representations that could confuse the algorithm.
    // xff\" and \xff\\ are multibyte sequences, so they don't get escaped. Other chars are escaped.
    string_view s = "This is \\ a string \xff\\ with a weird \xff\" encoding \"";
    std::string output = "abc";

    auto ec = escape_string(s, {test::ff_charset, true}, quoting_context::double_quote, output);

    BOOST_TEST(ec == error_code());
    BOOST_TEST(output == "This is \\\\ a string \xff\\ with a weird \xff\" encoding \\\"");
}

//
// Escaping doubling quotes
//
BOOST_AUTO_TEST_CASE(quotes_utf8mb4_valid)
{
    std::string all_ascii_noescape;
    all_ascii_noescape.reserve(0x80);
    for (int i = 0; i < 0x80; ++i)
    {
        char c = static_cast<char>(i);
        if (c != '"')
            all_ascii_noescape.push_back(c);
    }

    struct
    {
        string_view name;
        string_view input;
        string_view expected;
    } test_cases[] = {
        {"empty",              "",                                      ""                                       },
        {"no_escape_ascii",    "this is A test string",                 "this is A test string"                  },
        {"all_noescape_ascii", all_ascii_noescape,                      all_ascii_noescape                       },
        {"escape_quotes",      R"(this has "dquotes")",                 R"(this has ""dquotes"")"                },
        {"other_escape_chars", R"('squotes' "and" `backticks`\)",       R"('squotes' ""and"" `backticks`\)"      },
        {"all_escape_chars",   R"(""""")",                              R"("""""""""")"                          },
        {"start_escape_char",  "\"abc",                                 R"(""abc)"                               },
        {"end_escape_char",    "abc\"",                                 R"(abc"")"                               },
        {"single_escape_char", "\"",                                    R"("")"                                  },
        {"utf8_2byte",         "2byte \" \xc3\xb1 UTF-8\\ \xc3\xb2 \\", "2byte \"\" \xc3\xb1 UTF-8\\ \xc3\xb2 \\"},
        {"utf8_3byte",         "3byte \"\xef\xbf\xbf UTF-8\"",          "3byte \"\"\xef\xbf\xbf UTF-8\"\""       },
        {"utf8_4byte",         "4byte \"\xf0\x90\x80\x80 UTF-8\"",      "4byte \"\"\xf0\x90\x80\x80 UTF-8\"\""   },
 // Some typical injection payloads
        {"injection_1",        R"(\\)",                                 R"(\\)"                                  },
        {"injection_2",        R"(' or ")",                             R"(' or "")"                             },
        {"injection_3",        R"(-- or #)",                            R"(-- or #)"                             },
        {"injection_4",        R"(' OR '1)",                            R"(' OR '1)"                             },
        {"injection_5",        R"(' OR 1 -- -)",                        R"(' OR 1 -- -)"                         },
        {"injection_6",        R"(" OR "" = ")",                        R"("" OR """" = "")"                     },
        {"injection_7",        R"(" OR 1 = 1 -- -)",                    R"("" OR 1 = 1 -- -)"                    },
        {"injection_8",        R"(' OR '' = ')",                        R"(' OR '' = ')"                         },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            std::string output = "abc";
            auto
                ec = escape_string(tc.input, {utf8mb4_charset, false}, quoting_context::double_quote, output);
            BOOST_TEST(ec == error_code());
            BOOST_TEST(output == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(quotes_quoting_contexts)
{
    // Test through all quoting contexts
    string_view input = R"(A "string" that 'contains' some `quotes` \'"`)";

    struct
    {
        string_view name;
        quoting_context quot_ctx;
        string_view expected;
    } test_cases[] = {
        {"dquote",   quoting_context::double_quote, R"(A ""string"" that 'contains' some `quotes` \'""`)"},
        {"squote",   quoting_context::single_quote, R"(A "string" that ''contains'' some `quotes` \''"`)"},
        {"backtick", quoting_context::backtick,     R"(A "string" that 'contains' some ``quotes`` \'"``)"},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            std::string output = "abc";
            auto ec = escape_string(input, {utf8mb4_charset, false}, tc.quot_ctx, output);
            BOOST_TEST(ec == error_code());
            BOOST_TEST(output == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(quotes_utf8mb4_invalid)
{
    // \xc3\\ is an attempt to smuggle a backslash as an invalid 2 byte UTF8 sequence
    std::string output = "abc";
    auto ec = escape_string(
        "This \"has\" invalid \xc3\\ chars",
        {utf8mb4_charset, false},
        quoting_context::double_quote,
        output
    );
    BOOST_TEST(ec == client_errc::invalid_encoding);
}

BOOST_AUTO_TEST_CASE(quotes_multibyte_ascii_compatible_chars)
{
    // Edge cases for encodings allowing character representations that could confuse the algorithm.
    // \xff\" is a multibyte sequence, so it doesn't get escaped. Other chars are escaped.
    string_view s = "This is \" a string \xfe\" with a weird \xff\" encoding \"";
    std::string output = "abc";

    auto ec = escape_string(s, {test::ff_charset, false}, quoting_context::double_quote, output);

    BOOST_TEST(ec == error_code());
    BOOST_TEST(output == "This is \"\" a string \xfe\"\" with a weird \xff\" encoding \"\"");
}

BOOST_AUTO_TEST_CASE(parameter_coverage)
{
    // Test that the different combination of parameters dispatch
    // to the algorithm they should (backslashes or quotes)
    string_view input = "This \"has\" 'squotes'\n, `backticks`, and \\";

    struct
    {
        string_view name;
        bool backslash_escapes;
        quoting_context quot_ctx;
        string_view expected;
    } test_cases[] = {
        {"escapes_dquotes",
         true,  quoting_context::double_quote,
         "This \\\"has\\\" \\'squotes\\'\\n, `backticks`, and \\\\"},
        {"escapes_squotes",
         true,  quoting_context::single_quote,
         "This \\\"has\\\" \\'squotes\\'\\n, `backticks`, and \\\\"},
        {"escapes_backticks",
         true,  quoting_context::backtick,
         "This \"has\" 'squotes'\n, ``backticks``, and \\"         },

        {"no_escapes_dquotes",
         false, quoting_context::double_quote,
         "This \"\"has\"\" 'squotes'\n, `backticks`, and \\"       },
        {"no_escapes_squotes",
         false, quoting_context::single_quote,
         "This \"has\" ''squotes''\n, `backticks`, and \\"         },
        {"no_escapes_backticks",
         false, quoting_context::backtick,
         "This \"has\" 'squotes'\n, ``backticks``, and \\"         },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            std::string output = "abc";
            auto ec = escape_string(input, {utf8mb4_charset, tc.backslash_escapes}, tc.quot_ctx, output);
            BOOST_TEST(ec == error_code());
            BOOST_TEST(output == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(other_string_types)
{
    // Spotcheck: escape_string can be used with string types != std::string
    std::basic_string<char, std::char_traits<char>, test::custom_allocator<char>> output = "abc";
    auto ec = escape_string("some 'value'", {utf8mb4_charset, true}, quoting_context::single_quote, output);
    BOOST_TEST(ec == error_code());
    BOOST_TEST(output == R"(some \'value\')");
}

BOOST_AUTO_TEST_SUITE_END()
