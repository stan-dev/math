//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/character_set.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include <cassert>
#include <cstddef>

namespace mysql = boost::mysql;

namespace {

//[charsets_next_char
// next_char must interpret input as a string encoded according to the
// utf8mb4 character set and return the size of the first character,
// or 0 if the byte sequence does not represent a valid character.
// It must not throw exceptions.
std::size_t utf8mb4_next_char(boost::span<const unsigned char> input)
{
    // Input strings are never empty - they always have 1 byte, at least.
    assert(!input.empty());

    // In UTF8, we need to look at the first byte to know the character's length
    auto first_char = input[0];

    if (first_char < 0x80)
    {
        // 0x00 to 0x7F: ASCII range. The character is 1 byte long
        return 1;
    }
    else if (first_char <= 0xc1)
    {
        // 0x80 to 0xc1: invalid. No UTF8 character starts with such a byte
        return 0;
    }
    else if (first_char <= 0xdf)
    {
        // 0xc2 to 0xdf: two byte characters.
        // It's vital that we check that the characters are valid. Otherwise, vulnerabilities can arise.

        // Check that the string has enough bytes
        if (input.size() < 2u)
            return 0;

        // The second byte must be between 0x80 and 0xbf. Otherwise, the character is invalid
        // Do not skip this check - otherwise escaping will yield invalid results
        if (input[1] < 0x80 || input[1] > 0xbf)
            return 0;

        // Valid, 2 byte character
        return 2;
    }
    // Omitted: 3 and 4 byte long characters
    else
    {
        return 0;
    }
}
//]

BOOST_AUTO_TEST_CASE(section_charsets)
{
    {
        // Verify that utf8mb4_next_char can be used in a character_set
        mysql::character_set charset{"utf8mb4", utf8mb4_next_char};

        // It works for valid input
        unsigned char buff_valid[] = {0xc3, 0xb1, 0x50};
        BOOST_TEST(charset.next_char({buff_valid, sizeof(buff_valid)}) == 2u);

        // It works for invalid input
        unsigned char buff_invalid[] = {0xc3, 0xff, 0x50};
        BOOST_TEST(charset.next_char({buff_invalid, sizeof(buff_invalid)}) == 0u);
    }
}

}  // namespace