/*=============================================================================
    Copyright (c) 2020 Nikita Kniazev

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#include <boost/spirit/home/support/char_encoding/standard_wide.hpp>
#include <boost/spirit/home/support/char_encoding/unicode.hpp>

#include <boost/core/lightweight_test.hpp>

#if defined(_MSC_VER) && _MSC_VER < 1700
# pragma warning(disable: 4428) // universal-character-name encountered in source
#endif

int main()
{
    {
        using boost::spirit::char_encoding::standard_wide;
        BOOST_TEST_EQ(standard_wide::toucs4(L'\uFFFF'), 0xFFFFu);
        BOOST_TEST_EQ(standard_wide::toucs4(L'\u7FFF'), 0x7FFFu);
        BOOST_TEST_EQ(standard_wide::toucs4(L'\u0024'), 0x0024u);
    }

    {   // Unicode major categories
        using namespace boost::spirit::char_encoding;
        BOOST_TEST(unicode::is_letter(0x037Fu));
        BOOST_TEST(!unicode::is_letter(0x065Fu));
        BOOST_TEST(unicode::is_mark(0x065Fu));
        BOOST_TEST(!unicode::is_mark(0x0DE6u));
        BOOST_TEST(unicode::is_number(0x0DE6u));
        BOOST_TEST(!unicode::is_number(0x3000u));
        BOOST_TEST(unicode::is_separator(0x3000u));
        BOOST_TEST(!unicode::is_separator(0x0604u));
        BOOST_TEST(unicode::is_other(0x0604u));
        BOOST_TEST(!unicode::is_other(0x2E3Cu));
        BOOST_TEST(unicode::is_punctuation(0x2E3Cu));
        BOOST_TEST(!unicode::is_punctuation(0x26CEu));
        BOOST_TEST(unicode::is_symbol(0x26CEu));
        BOOST_TEST(!unicode::is_symbol(0x037Fu));
    }

    {   // Unicode general categories
        using namespace boost::spirit::char_encoding;

        BOOST_TEST(unicode::is_uppercase_letter(0x037Fu));
        BOOST_TEST(!unicode::is_uppercase_letter(0x065Fu));
        BOOST_TEST(unicode::is_lowercase_letter(0xABBAu));
        BOOST_TEST(!unicode::is_lowercase_letter(0x1ABEu));
        BOOST_TEST(unicode::is_titlecase_letter(0x1F88u));
        BOOST_TEST(!unicode::is_titlecase_letter(0x1BE7u));
        BOOST_TEST(unicode::is_modifier_letter(0x08C9u));
        BOOST_TEST(!unicode::is_modifier_letter(0xA9F0u));
        BOOST_TEST(unicode::is_other_letter(0x08C8u));
        BOOST_TEST(!unicode::is_other_letter(0x12463u));

        BOOST_TEST(unicode::is_nonspacing_mark(0x065Fu));
        BOOST_TEST(!unicode::is_nonspacing_mark(0x16B5Bu));
        BOOST_TEST(unicode::is_enclosing_mark(0x1ABEu));
        BOOST_TEST(!unicode::is_enclosing_mark(0x06DEu));
        BOOST_TEST(unicode::is_spacing_mark(0x1BE7u));
        BOOST_TEST(!unicode::is_spacing_mark(0x2000u));

        BOOST_TEST(unicode::is_decimal_number(0xA9F0u));
        BOOST_TEST(!unicode::is_decimal_number(0x2028u));
        BOOST_TEST(unicode::is_letter_number(0x12463u));
        BOOST_TEST(!unicode::is_letter_number(0x2029u));
        BOOST_TEST(unicode::is_other_number(0x16B5Bu));
        BOOST_TEST(!unicode::is_other_number(0x0604u));

        BOOST_TEST(unicode::is_space_separator(0x2000u));
        BOOST_TEST(!unicode::is_space_separator(0x10FFFDu));
        BOOST_TEST(unicode::is_line_separator(0x2028u));
        BOOST_TEST(!unicode::is_line_separator(0xDDDDu));
        BOOST_TEST(unicode::is_paragraph_separator(0x2029u));
        BOOST_TEST(!unicode::is_paragraph_separator(0x1FFEu));

        BOOST_TEST(unicode::is_control(0x0000u));
        BOOST_TEST(!unicode::is_control(0x037Fu));
        BOOST_TEST(unicode::is_format(0x0604u));
        BOOST_TEST(!unicode::is_format(0x0606u));
        BOOST_TEST(unicode::is_private_use(0x10FFFDu));
        BOOST_TEST(!unicode::is_private_use(0xABBAu));
        BOOST_TEST(unicode::is_surrogate(0xDDDDu));
        BOOST_TEST(!unicode::is_surrogate(0x1F88u));
        BOOST_TEST(unicode::is_unassigned(0x1FFFu));
        BOOST_TEST(!unicode::is_unassigned(0x1FFEu));

        BOOST_TEST(unicode::is_dash_punctuation(0x2E3Au));
        BOOST_TEST(!unicode::is_dash_punctuation(0x08C9u));
        BOOST_TEST(unicode::is_open_punctuation(0x2E42u));
        BOOST_TEST(!unicode::is_open_punctuation(0x08C8u));
        BOOST_TEST(unicode::is_close_punctuation(0x2E56u));
        BOOST_TEST(!unicode::is_close_punctuation(0x065Fu));
        BOOST_TEST(unicode::is_connector_punctuation(0x203Fu));
        BOOST_TEST(!unicode::is_connector_punctuation(0x1ABEu));
        BOOST_TEST(unicode::is_other_punctuation(0x00B6u));
        BOOST_TEST(!unicode::is_other_punctuation(0x1BE7u));
        BOOST_TEST(unicode::is_initial_punctuation(0x00ABu));
        BOOST_TEST(!unicode::is_initial_punctuation(0xA9F0u));
        BOOST_TEST(unicode::is_final_punctuation(0x00BBu));
        BOOST_TEST(!unicode::is_final_punctuation(0x12463u));

        BOOST_TEST(unicode::is_math_symbol(0x27CBu));
        BOOST_TEST(!unicode::is_math_symbol(0x16B5Bu));
        BOOST_TEST(unicode::is_currency_symbol(0x11FDDu));
        BOOST_TEST(!unicode::is_currency_symbol(0x2000u));
        BOOST_TEST(unicode::is_modifier_symbol(0xAB5Bu));
        BOOST_TEST(!unicode::is_modifier_symbol(0x2028u));
        BOOST_TEST(unicode::is_other_symbol(0x27BFu));
        BOOST_TEST(!unicode::is_other_symbol(0x2029u));
    }

    {   // Unicode derived categories
        using namespace boost::spirit::char_encoding;
        BOOST_TEST(unicode::is_alphabetic(0x0555u));
        BOOST_TEST(!unicode::is_alphabetic(0x0557u));
        BOOST_TEST(unicode::is_uppercase(0x10410u));
        BOOST_TEST(!unicode::is_uppercase(0x10430u));
        BOOST_TEST(unicode::is_lowercase(0x00AAu));
        BOOST_TEST(!unicode::is_lowercase(0x00ABu));
        BOOST_TEST(unicode::is_white_space(0x2002u));
        BOOST_TEST(!unicode::is_white_space(0x200Bu));
        BOOST_TEST(unicode::is_hex_digit(0xFF26u));
        BOOST_TEST(!unicode::is_hex_digit(0xFF27u));
        BOOST_TEST(unicode::is_noncharacter_code_point(0x10FFFEu));
        BOOST_TEST(!unicode::is_noncharacter_code_point(0x10FFFDu));
        BOOST_TEST(unicode::is_default_ignorable_code_point(0xE0FFFu));
        BOOST_TEST(!unicode::is_default_ignorable_code_point(0xE1000u));
    }

    {   // Unicode scripts
        using namespace boost::spirit::char_encoding;
        BOOST_TEST(unicode::is_arabic(0x060Du));
        BOOST_TEST(!unicode::is_arabic(0xE000u));
        BOOST_TEST(unicode::is_braille(0x2828u));
        BOOST_TEST(!unicode::is_braille(0x2728u));
        BOOST_TEST(unicode::is_toto(0x1E290u));
        BOOST_TEST(!unicode::is_toto(0x1E2AFu));
        BOOST_TEST(unicode::is_inherited(0x0300u));
        BOOST_TEST(!unicode::is_inherited(0x02FFu));
        BOOST_TEST(unicode::is_common(0xE0001u));
        BOOST_TEST(!unicode::is_common(0xE0000u));
        BOOST_TEST(unicode::is_unknown(0xA63Fu));
        BOOST_TEST(unicode::is_unknown(0xD800u));
        BOOST_TEST(unicode::is_unknown(0xE000u));
        BOOST_TEST(unicode::is_unknown(0x10FFFFu));
        BOOST_TEST(!unicode::is_unknown(0xE0001u));
    }

    return boost::report_errors();
}
