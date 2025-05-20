//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/impl/internal/call_next_char.hpp>

#include <boost/test/tools/context.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <string>

#include "test_common/create_basic.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_AUTO_TEST_SUITE(test_character_set)

// Helper
static std::size_t call_next_char(const character_set& charset, string_view s) noexcept
{
    return detail::call_next_char(charset, s.data(), s.data() + s.size());
}

BOOST_AUTO_TEST_CASE(utf8mb4_single_byte_valid)
{
    for (int i = 0; i < 0x80; ++i)
    {
        BOOST_TEST_CONTEXT(i)
        {
            // Exactly the required space
            char str[2]{static_cast<char>(i), '\0'};
            auto actual_len = detail::call_next_char(utf8mb4_charset, str, str + 1);
            BOOST_TEST(actual_len == 1u);

            // Extra space
            actual_len = detail::call_next_char(utf8mb4_charset, str, str + 2);
            BOOST_TEST(actual_len == 1u);
        }
    }
}

BOOST_AUTO_TEST_CASE(utf8mb4_multibyte_valid)
{
    // Cases are structured depending on the first byte.
    // We cover all possible first bytes, with some casuistic for each value.
    struct
    {
        string_view name;
        string_view input;
        std::size_t expected;
    } test_cases[] = {
  // 2 byte characters. We perform some extra tests for c2 and c3
        {"c2 min (U+0080)",          "\xc2\x80",         2u},
        {"c2 reg (U+0095)",          "\xc2\x95",         2u},
        {"c2 max (U+00BF)",          "\xc2\xbf",         2u},
        {"c3 min (U+00C0)",          "\xc3\x80",         2u},
        {"c3 reg (U+00E7)",          "\xc3\xa7",         2u},
        {"c3 max (U+00FF)",          "\xc3\xbf",         2u},

 // c4-df behave the same as c2-c3, so only min and max
        {"c4 min (U+0100)",          "\xc4\x80",         2u},
        {"c4 max (U+013F)",          "\xc4\xbf",         2u},
        {"c5 min (U+0140)",          "\xc5\x80",         2u},
        {"c5 max (U+017F)",          "\xc5\xbf",         2u},
        {"c6 min (U+0180)",          "\xc6\x80",         2u},
        {"c6 max (U+01BF)",          "\xc6\xbf",         2u},
        {"c7 min (U+01C0)",          "\xc7\x80",         2u},
        {"c7 max (U+01FF)",          "\xc7\xbf",         2u},
        {"c8 min (U+0200)",          "\xc8\x80",         2u},
        {"c8 max (U+023F)",          "\xc8\xbf",         2u},
        {"c9 min (U+0240)",          "\xc9\x80",         2u},
        {"c9 max (U+027F)",          "\xc9\xbf",         2u},
        {"ca min (U+0280)",          "\xca\x80",         2u},
        {"ca max (U+02BF)",          "\xca\xbf",         2u},
        {"cb min (U+02C0)",          "\xcb\x80",         2u},
        {"cb max (U+02FF)",          "\xcb\xbf",         2u},
        {"cc min (U+0300)",          "\xcc\x80",         2u},
        {"cc max (U+033F)",          "\xcc\xbf",         2u},
        {"cd min (U+0340)",          "\xcd\x80",         2u},
        {"cd max (U+037F)",          "\xcd\xbf",         2u},
        {"ce min (U+0380)",          "\xce\x80",         2u},
        {"ce max (U+03BF)",          "\xce\xbf",         2u},
        {"cf min (U+03C0)",          "\xcf\x80",         2u},
        {"cf max (U+03FF)",          "\xcf\xbf",         2u},
        {"d0 min (U+0400)",          "\xd0\x80",         2u},
        {"d0 max (U+043F)",          "\xd0\xbf",         2u},
        {"d1 min (U+0440)",          "\xd1\x80",         2u},
        {"d1 max (U+047F)",          "\xd1\xbf",         2u},
        {"d2 min (U+0480)",          "\xd2\x80",         2u},
        {"d2 max (U+04BF)",          "\xd2\xbf",         2u},
        {"d3 min (U+04C0)",          "\xd3\x80",         2u},
        {"d3 max (U+04FF)",          "\xd3\xbf",         2u},
        {"d4 min (U+0500)",          "\xd4\x80",         2u},
        {"d4 max (U+053F)",          "\xd4\xbf",         2u},
        {"d5 min (U+0540)",          "\xd5\x80",         2u},
        {"d5 max (U+057F)",          "\xd5\xbf",         2u},
        {"d6 min (U+0580)",          "\xd6\x80",         2u},
        {"d6 max (U+05BF)",          "\xd6\xbf",         2u},
        {"d7 min (U+05C0)",          "\xd7\x80",         2u},
        {"d7 max (U+05FF)",          "\xd7\xbf",         2u},
        {"d8 min (U+0600)",          "\xd8\x80",         2u},
        {"d8 max (U+063F)",          "\xd8\xbf",         2u},
        {"d9 min (U+0640)",          "\xd9\x80",         2u},
        {"d9 max (U+067F)",          "\xd9\xbf",         2u},
        {"da min (U+0680)",          "\xda\x80",         2u},
        {"da max (U+06BF)",          "\xda\xbf",         2u},
        {"db min (U+06C0)",          "\xdb\x80",         2u},
        {"db max (U+06FF)",          "\xdb\xbf",         2u},
        {"dc min (U+0700)",          "\xdc\x80",         2u},
        {"dc max (U+073F)",          "\xdc\xbf",         2u},
        {"dd min (U+0740)",          "\xdd\x80",         2u},
        {"dd max (U+077F)",          "\xdd\xbf",         2u},
        {"de min (U+0780)",          "\xde\x80",         2u},
        {"de max (U+07BF)",          "\xde\xbf",         2u},
        {"df min (U+07C0)",          "\xdf\x80",         2u},
        {"df max (U+07FF)",          "\xdf\xbf",         2u},

 // 3 byte characters. We perform some extra tests for e0 and e1
        {"e0 min min (U+0800)",      "\xe0\xa0\x80",     3u},
        {"e0 min reg (U+0835)",      "\xe0\xa0\xb5",     3u},
        {"e0 min max (U+083F)",      "\xe0\xa0\xbf",     3u},
        {"e0 reg min (U+0900)",      "\xe0\xa4\x80",     3u},
        {"e0 reg reg (U+0920)",      "\xe0\xa4\xa0",     3u},
        {"e0 reg max (U+093F)",      "\xe0\xa4\xbf",     3u},
        {"e0 max min (U+0FC0)",      "\xe0\xbf\x80",     3u},
        {"e0 max reg (U+0FE0)",      "\xe0\xbf\xa0",     3u},
        {"e0 max max (U+0FFF)",      "\xe0\xbf\xbf",     3u},

        {"e1 min min (U+1000)",      "\xe1\x80\x80",     3u},
        {"e1 min max (U+103F)",      "\xe1\x80\xbf",     3u},
        {"e1 max min (U+1FC0)",      "\xe1\xbf\x80",     3u},
        {"e1 max max (U+1FFF)",      "\xe1\xbf\xbf",     3u},

 // e2-ec behave like e1
        {"e2 min (U+2000)",          "\xe2\x80\x80",     3u},
        {"e2 max (U+2FFF)",          "\xe2\xbf\xbf",     3u},
        {"e3 min (U+3000)",          "\xe3\x80\x80",     3u},
        {"e3 max (U+3FFF)",          "\xe3\xbf\xbf",     3u},
        {"e4 min (U+4000)",          "\xe4\x80\x80",     3u},
        {"e4 max (U+4FFF)",          "\xe4\xbf\xbf",     3u},
        {"e5 min (U+5000)",          "\xe5\x80\x80",     3u},
        {"e5 max (U+5FFF)",          "\xe5\xbf\xbf",     3u},
        {"e6 min (U+6000)",          "\xe6\x80\x80",     3u},
        {"e6 max (U+6FFF)",          "\xe6\xbf\xbf",     3u},
        {"e7 min (U+7000)",          "\xe7\x80\x80",     3u},
        {"e7 max (U+7FFF)",          "\xe7\xbf\xbf",     3u},
        {"e8 min (U+8000)",          "\xe8\x80\x80",     3u},
        {"e8 max (U+8FFF)",          "\xe8\xbf\xbf",     3u},
        {"e9 min (U+9000)",          "\xe9\x80\x80",     3u},
        {"e9 max (U+9FFF)",          "\xe9\xbf\xbf",     3u},
        {"ea min (U+A000)",          "\xea\x80\x80",     3u},
        {"ea max (U+AFFF)",          "\xea\xbf\xbf",     3u},
        {"eb min (U+B000)",          "\xeb\x80\x80",     3u},
        {"eb max (U+BFFF)",          "\xeb\xbf\xbf",     3u},
        {"ec min (U+C000)",          "\xec\x80\x80",     3u},
        {"ec max (U+CFFF)",          "\xec\xbf\xbf",     3u},

 // ed is different because of surrogates (code points U+D800 to U+DFFF)
        {"ed min (U+D000)",          "\xed\x80\x80",     3u},
        {"ed reg (U+D631)",          "\xed\x98\xb1",     3u},
        {"ed max (U+D7FF)",          "\xed\x9f\xbf",     3u},

 // ee-ef behave like e1
        {"ee min (U+E000)",          "\xee\x80\x80",     3u},
        {"ee max (U+EFFF)",          "\xee\xbf\xbf",     3u},
        {"ef min (U+F000)",          "\xef\x80\x80",     3u},
        {"ef max (U+FFFF)",          "\xef\xbf\xbf",     3u},

 // 4 byte characters - we perform some extra testing for f0
        {"f0 min min min (U+10000)", "\xf0\x90\x80\x80", 4u},
        {"f0 min min max (U+1003F)", "\xf0\x90\x80\xbf", 4u},
        {"f0 min max min (U+10FC0)", "\xf0\x90\xbf\x80", 4u},
        {"f0 min max max (U+10FFF)", "\xf0\x90\xbf\xbf", 4u},
        {"f0 max min min (U+3F000)", "\xf0\xbf\x80\x80", 4u},
        {"f0 max min max (U+3F03F)", "\xf0\xbf\x80\xbf", 4u},
        {"f0 max max min (U+3FFC0)", "\xf0\xbf\xbf\x80", 4u},
        {"f0 max max max (U+3FFFF)", "\xf0\xbf\xbf\xbf", 4u},

        {"f1 min (U+40000)",         "\xf1\x80\x80\x80", 4u},
        {"f1 max (U+7FFFF)",         "\xf1\xbf\xbf\xbf", 4u},
        {"f2 min (U+80000)",         "\xf2\x80\x80\x80", 4u},
        {"f2 max (U+BFFFF)",         "\xf2\xbf\xbf\xbf", 4u},
        {"f3 min (U+C0000)",         "\xf3\x80\x80\x80", 4u},
        {"f3 max (U+FFFFF)",         "\xf3\xbf\xbf\xbf", 4u},

 // The last allowable code point is U+10FFFF
        {"f4 min (U+100000)",        "\xf4\x80\x80\x80", 4u},
        {"f4 max (U+10FFFF)",        "\xf4\x8f\xbf\xbf", 4u},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Exactly the required space
            auto actual_len = call_next_char(utf8mb4_charset, tc.input);
            BOOST_TEST(actual_len == tc.expected);

            // Extra space
            auto extra_space_input = std::string(tc.input) + "abc";
            actual_len = call_next_char(utf8mb4_charset, extra_space_input);
            BOOST_TEST(actual_len == tc.expected);

            // Not enough space (end of data before the end of the byte sequence)
            auto not_enough_input = tc.input.substr(1);
            actual_len = call_next_char(utf8mb4_charset, not_enough_input);
            BOOST_TEST(actual_len == 0u);
        }
    }
}

BOOST_AUTO_TEST_CASE(utf8mb4_invalid_start_byte)
{
    constexpr unsigned char bytes[] = {
        0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
        0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
        0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
        0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
        0xc0, 0xc1, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff,
    };

    for (const auto b : bytes)
    {
        BOOST_TEST_CONTEXT(+b)
        {
            auto input = static_cast<char>(b);
            auto size = detail::call_next_char(utf8mb4_charset, &input, &input + 1);
            BOOST_TEST(size == 0u);
        }
    }
}

BOOST_AUTO_TEST_CASE(utf8mb4_invalid_continuation)
{
    struct
    {
        string_view name;
        string_view input;
    } test_cases[] = {
  // 2 byte characters
        {"c2 zero",             makesv("\xc2\x00")        },
        {"c2 ltmin",            "\xc2\x7f"                },
        {"c2 gtmax",            "\xc2\xc0"                },
        {"c2 max",              "\xc2\xff"                },
        {"c3 ltmin",            "\xc3\x7f"                },
        {"c3 gtmax",            "\xc3\xc0"                },
        {"c4 ltmin",            "\xc4\x7f"                },
        {"c4 gtmax",            "\xc4\xc0"                },
        {"c5 ltmin",            "\xc5\x7f"                },
        {"c5 gtmax",            "\xc5\xc0"                },
        {"c6 ltmin",            "\xc6\x7f"                },
        {"c6 gtmax",            "\xc6\xc0"                },
        {"c7 ltmin",            "\xc7\x7f"                },
        {"c7 gtmax",            "\xc7\xc0"                },
        {"c8 ltmin",            "\xc8\x7f"                },
        {"c8 gtmax",            "\xc8\xc0"                },
        {"c9 ltmin",            "\xc9\x7f"                },
        {"c9 gtmax",            "\xc9\xc0"                },
        {"ca ltmin",            "\xca\x7f"                },
        {"ca gtmax",            "\xca\xc0"                },
        {"cb ltmin",            "\xcb\x7f"                },
        {"cb gtmax",            "\xcb\xc0"                },
        {"cc ltmin",            "\xcc\x7f"                },
        {"cc gtmax",            "\xcc\xc0"                },
        {"cd ltmin",            "\xcd\x7f"                },
        {"cd gtmax",            "\xcd\xc0"                },
        {"ce ltmin",            "\xce\x7f"                },
        {"ce gtmax",            "\xce\xc0"                },
        {"cf ltmin",            "\xcf\x7f"                },
        {"cf gtmax",            "\xcf\xc0"                },
        {"d0 ltmin",            "\xd0\x7f"                },
        {"d0 gtmax",            "\xd0\xc0"                },
        {"d1 ltmin",            "\xd1\x7f"                },
        {"d1 gtmax",            "\xd1\xc0"                },
        {"d2 ltmin",            "\xd2\x7f"                },
        {"d2 gtmax",            "\xd2\xc0"                },
        {"d3 ltmin",            "\xd3\x7f"                },
        {"d3 gtmax",            "\xd3\xc0"                },
        {"d4 ltmin",            "\xd4\x7f"                },
        {"d4 gtmax",            "\xd4\xc0"                },
        {"d5 ltmin",            "\xd5\x7f"                },
        {"d5 gtmax",            "\xd5\xc0"                },
        {"d6 ltmin",            "\xd6\x7f"                },
        {"d6 gtmax",            "\xd6\xc0"                },
        {"d7 ltmin",            "\xd7\x7f"                },
        {"d7 gtmax",            "\xd7\xc0"                },
        {"d8 ltmin",            "\xd8\x7f"                },
        {"d8 gtmax",            "\xd8\xc0"                },
        {"d9 ltmin",            "\xd9\x7f"                },
        {"d9 gtmax",            "\xd9\xc0"                },
        {"da ltmin",            "\xda\x7f"                },
        {"da gtmax",            "\xda\xc0"                },
        {"db ltmin",            "\xdb\x7f"                },
        {"db gtmax",            "\xdb\xc0"                },
        {"dc ltmin",            "\xdc\x7f"                },
        {"dc gtmax",            "\xdc\xc0"                },
        {"dd ltmin",            "\xdd\x7f"                },
        {"dd gtmax",            "\xdd\xc0"                },
        {"de ltmin",            "\xde\x7f"                },
        {"de gtmax",            "\xde\xc0"                },
        {"df ltmin",            "\xdf\x7f"                },
        {"df gtmax",            "\xdf\xc0"                },

 // 3 byte chars (e0 is special)
        {"e0 ltmin ok",         "\xe0\x9f\x91"            },
        {"e0 gtmax ok",         "\xe0\xc0\x91"            },
        {"e0 ok ltmin",         "\xe0\xa0\x7F"            },
        {"e0 ok gtmax",         "\xe0\xa0\xc0"            },

        {"e1 ltmin ok",         "\xe1\x7f\x91"            },
        {"e1 gtmax ok",         "\xe1\xc0\x91"            },
        {"e1 ok ltmin",         "\xe1\xa0\x7F"            },
        {"e1 ok gtmax",         "\xe1\xa0\xc0"            },

        {"e2 ltmin ok",         "\xe2\x7f\x91"            },
        {"e2 gtmax ok",         "\xe2\xc0\x91"            },
        {"e2 ok ltmin",         "\xe2\xa0\x7F"            },
        {"e2 ok gtmax",         "\xe2\xa0\xc0"            },

        {"e3 ltmin ok",         "\xe3\x7f\x91"            },
        {"e3 gtmax ok",         "\xe3\xc0\x91"            },
        {"e3 ok ltmin",         "\xe3\xa0\x7F"            },
        {"e3 ok gtmax",         "\xe3\xa0\xc0"            },

        {"e4 ltmin ok",         "\xe4\x7f\x91"            },
        {"e4 gtmax ok",         "\xe4\xc0\x91"            },
        {"e4 ok ltmin",         "\xe4\xa0\x7F"            },
        {"e4 ok gtmax",         "\xe4\xa0\xc0"            },

        {"e5 ltmin ok",         "\xe5\x7f\x91"            },
        {"e5 gtmax ok",         "\xe5\xc0\x91"            },
        {"e5 ok ltmin",         "\xe5\xa0\x7F"            },
        {"e5 ok gtmax",         "\xe5\xa0\xc0"            },

        {"e6 ltmin ok",         "\xe6\x7f\x91"            },
        {"e6 gtmax ok",         "\xe6\xc0\x91"            },
        {"e6 ok ltmin",         "\xe6\xa0\x7F"            },
        {"e6 ok gtmax",         "\xe6\xa0\xc0"            },

        {"e7 ltmin ok",         "\xe7\x7f\x91"            },
        {"e7 gtmax ok",         "\xe7\xc0\x91"            },
        {"e7 ok ltmin",         "\xe7\xa0\x7F"            },
        {"e7 ok gtmax",         "\xe7\xa0\xc0"            },

        {"e8 ltmin ok",         "\xe8\x7f\x91"            },
        {"e8 gtmax ok",         "\xe8\xc0\x91"            },
        {"e8 ok ltmin",         "\xe8\xa0\x7F"            },
        {"e8 ok gtmax",         "\xe8\xa0\xc0"            },

        {"e9 ltmin ok",         "\xe9\x7f\x91"            },
        {"e9 gtmax ok",         "\xe9\xc0\x91"            },
        {"e9 ok ltmin",         "\xe9\xa0\x7F"            },
        {"e9 ok gtmax",         "\xe9\xa0\xc0"            },

        {"ea ltmin ok",         "\xea\x7f\x91"            },
        {"ea gtmax ok",         "\xea\xc0\x91"            },
        {"ea ok ltmin",         "\xea\xa0\x7F"            },
        {"ea ok gtmax",         "\xea\xa0\xc0"            },

        {"eb ltmin ok",         "\xeb\x7f\x91"            },
        {"eb gtmax ok",         "\xeb\xc0\x91"            },
        {"eb ok ltmin",         "\xeb\xa0\x7F"            },
        {"eb ok gtmax",         "\xeb\xa0\xc0"            },

        {"ec ltmin ok",         "\xec\x7f\x91"            },
        {"ec gtmax ok",         "\xec\xc0\x91"            },
        {"ec ok ltmin",         "\xec\xa0\x7F"            },
        {"ec ok gtmax",         "\xec\xa0\xc0"            },

 // ed is special because it includes surrogates
        {"ed ltmin ok",         "\xed\x7f\x91"            },
        {"ed gtmax ok",         "\xed\xc0\x91"            },
        {"ed ok ltmin",         "\xed\xa0\x7F"            },
        {"ed ok gtmax",         "\xed\xa0\xc0"            },
        {"ed surrogate min",    "\xed\xa0\x80"            },
        {"ed surrogate reg",    "\xed\xa1\x92"            },
        {"ed surrogate max",    "\xed\xbf\xbf"            },

 // ee and ef behave like e1
        {"ee ltmin ok",         "\xee\x7f\x91"            },
        {"ee gtmax ok",         "\xee\xc0\x91"            },
        {"ee ok ltmin",         "\xee\xa0\x7F"            },
        {"ee ok gtmax",         "\xee\xa0\xc0"            },

        {"ef ltmin ok",         "\xef\x7f\x91"            },
        {"ef gtmax ok",         "\xef\xc0\x91"            },
        {"ef ok ltmin",         "\xef\xa0\x7F"            },
        {"ef ok gtmax",         "\xef\xa0\xc0"            },

 // 4 byte characters. f0 is special
        {"f0 ltmin ok ok",      "\xf0\x8f\x80\x80"        },
        {"f0 gtmax ok ok",      "\xf0\xc0\x80\x80"        },
        {"f0 ok ltmin ok",      "\xf0\xa1\x7f\xa3"        },
        {"f0 ok gtmax ok",      "\xf0\xa1\xc0\xa3"        },
        {"f0 ok ok ltmin",      "\xf0\xa1\xa2\x7f"        },
        {"f0 ok ok gtmax",      "\xf0\xa1\xa2\xc0"        },

        {"f1 ltmin ok ok",      "\xf1\x7f\x80\x80"        },
        {"f1 gtmax ok ok",      "\xf1\xc0\x80\x80"        },
        {"f1 ok ltmin ok",      "\xf1\xa1\x7f\xa3"        },
        {"f1 ok gtmax ok",      "\xf1\xa1\xc0\xa3"        },
        {"f1 ok ok ltmin",      "\xf1\xa1\xa2\x7f"        },
        {"f1 ok ok gtmax",      "\xf1\xa1\xa2\xc0"        },

        {"f2 ltmin ok ok",      "\xf2\x7f\x80\x80"        },
        {"f2 gtmax ok ok",      "\xf2\xc0\x80\x80"        },
        {"f2 ok ltmin ok",      "\xf2\xa1\x7f\xa3"        },
        {"f2 ok gtmax ok",      "\xf2\xa1\xc0\xa3"        },
        {"f2 ok ok ltmin",      "\xf2\xa1\xa2\x7f"        },
        {"f2 ok ok gtmax",      "\xf2\xa1\xa2\xc0"        },

        {"f3 ltmin ok ok",      "\xf3\x7f\x80\x80"        },
        {"f3 gtmax ok ok",      "\xf3\xc0\x80\x80"        },
        {"f3 ok ltmin ok",      "\xf3\xa1\x7f\xa3"        },
        {"f3 ok gtmax ok",      "\xf3\xa1\xc0\xa3"        },
        {"f3 ok ok ltmin",      "\xf3\xa1\xa2\x7f"        },
        {"f3 ok ok gtmax",      "\xf3\xa1\xa2\xc0"        },

 // f4 is also special because it's the end of the unicode range
        {"f4 ltmin ok ok",      "\xf4\x7f\x80\x80"        },
        {"f4 gtmax ok ok",      "\xf4\x90\x80\x80"        },
        {"f4 ok ltmin ok",      "\xf4\xa1\x7f\xa3"        },
        {"f4 ok gtmax ok",      "\xf4\xa1\xc0\xa3"        },
        {"f4 ok ok ltmin",      "\xf4\xa1\xa2\x7f"        },
        {"f4 ok ok gtmax",      "\xf4\xa1\xa2\xc0"        },

 // overlong characters
        {"overlong / 2byte",    "\xc0\xaf"                },
        {"overlong / 3byte",    "\xe0\x80\xaf"            },
        {"overlong / 4byte",    "\xf0\x80\x80\xaf"        },
        {"overlong / 5byte",    "\xf8\x80\x80\x80\xaf"    },
        {"overlong / 6byte",    "\xf8\x80\x80\x80\x80\xaf"},
        {"overlong U+007F",     "\xc1\xbf"                },
        {"overlong U+07FF",     "\xe0\x9f\xbf"            },
        {"overlong U+FFFF",     "\xf0\x8f\xbf\xbf"        },
        {"overlong U+001FFFFF", "\xf8\x87\xbf\xbf\xbf"    },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // add some extra continuation bytes, so we never fail because of lack of space
            auto input = std::string(tc.input) + "\x91\x91";
            auto size = call_next_char(utf8mb4_charset, input);
            BOOST_TEST(size == 0u);
        }
    }
}

BOOST_AUTO_TEST_CASE(ascii)
{
    // valid
    for (std::size_t i = 0; i <= 0x7f; ++i)
    {
        BOOST_TEST_CONTEXT(i)
        {
            char str[2]{static_cast<char>(i), '\0'};
            auto size = detail::call_next_char(ascii_charset, str, str + 2);
            BOOST_TEST(size == 1u);
        }
    }

    // invalid
    for (std::size_t i = 0x80; i <= 0xff; ++i)
    {
        BOOST_TEST_CONTEXT(i)
        {
            char str[2]{static_cast<char>(i), '\0'};
            auto size = detail::call_next_char(ascii_charset, str, str + 2);
            BOOST_TEST(size == 0u);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()