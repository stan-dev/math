/**
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/parser.hpp>
#include <boost/parser/transcode_view.hpp>

#include <boost/core/lightweight_test.hpp>


int main()
{

// doc_example)
{
    namespace bp = boost::parser;
    auto const street_parser = bp::string(u8"Tobias Straße");
    assert(!bp::parse("Tobias Strasse" | bp::as_utf32, street_parser));             // No match.
    assert(bp::parse("Tobias Strasse" | bp::as_utf32, bp::no_case[street_parser])); // Match!

    auto const alpha_parser = bp::no_case[bp::char_('a', 'z')];
    assert(bp::parse("a" | bp::as_utf32, bp::no_case[alpha_parser])); // Match!
    assert(bp::parse("B" | bp::as_utf32, bp::no_case[alpha_parser])); // Match!

    (void)street_parser;
    (void)alpha_parser;
}


using namespace boost::parser;

constexpr auto capital_sharp_s = u8"ẞ"; // U+1E9E
constexpr auto small_sharp_s = u8"ß";   // U+00DF
constexpr auto double_s = u8"sS";       // U+0073 U+0073

// basic)
{
    constexpr auto char_p = no_case[char_('a') | char_('B')];

    {
        auto const result = parse("a", char_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        auto const result = parse("A", char_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'A');
    }

    {
        auto const result = parse("b", char_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'b');
    }
    {
        auto const result = parse("B", char_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'B');
    }

    constexpr auto str_p = string("sOmE TeXT");

    {
        auto const result = parse("some text", str_p);
        BOOST_TEST(!result);
    }
    {
        auto const result = parse("some text", no_case[str_p]);
        BOOST_TEST(result);
        BOOST_TEST(*result == "some text");
    }
    {
        auto const result = parse("SomE tEXt", str_p);
        BOOST_TEST(!result);
    }
    {
        auto const result = parse("SomE tEXt", no_case[str_p]);
        BOOST_TEST(result);
        BOOST_TEST(*result == "SomE tEXt");
    }
}

// char_range)
{
    constexpr auto lower_alpha_p = char_('a', 'b');
    {
        BOOST_TEST(!parse("A", lower_alpha_p));
        auto const result = parse("a", lower_alpha_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }
    {
        auto const result = parse("A", no_case[lower_alpha_p]);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'A');
    }

    constexpr auto upper_alpha_p = char_('A', 'B');
    {
        BOOST_TEST(!parse("a", upper_alpha_p));
        auto const result = parse("A", upper_alpha_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'A');
    }
    {
        auto const result = parse("a", no_case[upper_alpha_p]);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'a');
    }

    constexpr auto weird_ss_p = char_(U'ß', U'ß');
    {
        BOOST_TEST(!parse(U"s", weird_ss_p));
        BOOST_TEST(!!parse(U"ß", weird_ss_p));
    }
}

// match_any_within_string)
{
    constexpr auto _trasse_p = no_case[char_(u8"_traße")];
    {
        auto const result = parse(U"ß", _trasse_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == U'ß');
    }
    {
        auto const result = parse(U"s", _trasse_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == U's');
    }
    {
        // Non-Unicode parsing fails to match, since 'ß' is not treated as a
        // single character.
        auto const result = parse("s", _trasse_p);
        BOOST_TEST(!result);
    }
    {
        auto const result = parse(U"S", _trasse_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == U'S');
    }
    {
        auto const result = parse("S", _trasse_p);
        BOOST_TEST(!result);
    }
    {
        auto const result = parse(U"t", _trasse_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == U't');
    }
    {
        auto const result = parse("t", _trasse_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 't');
    }
    {
        auto const result = parse(U"T", _trasse_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == U'T');
    }
    {
        auto const result = parse("T", _trasse_p);
        BOOST_TEST(result);
        BOOST_TEST(*result == 'T');
    }
    {
        auto const result = parse("X", _trasse_p);
        BOOST_TEST(!result);
    }
}

// symbol_table)
{
    // without mutation
    {
        symbols<int> const roman_numerals = {
            {"I", 1}, {"V", 5}, {"X", 10}, {"L", 50}, {"C", 100}};
        symbols<std::string> const named_strings = {
            {"I", "1"}, {"V", "5"}, {"X", "10"}, {"L", "50"}, {"C", "100"}};

        {
            auto const result = parse("I", no_case[roman_numerals]);
            BOOST_TEST(result);
            BOOST_TEST(*result == 1);
        }
        {
            auto const result = parse("i", no_case[roman_numerals]);
            BOOST_TEST(result);
            BOOST_TEST(*result == 1);
        }
        {
            auto const result = parse("I", no_case[named_strings]);
            BOOST_TEST(result);
            BOOST_TEST(*result == "1");
        }
        {
            auto const result = parse("i", no_case[named_strings]);
            BOOST_TEST(result);
            BOOST_TEST(*result == "1");
        }

        {
            auto const result = parse("L", no_case[roman_numerals]);
            BOOST_TEST(result);
            BOOST_TEST(*result == 50);
        }
        {
            auto const result = parse("l", no_case[roman_numerals]);
            BOOST_TEST(result);
            BOOST_TEST(*result == 50);
        }
        {
            auto const result = parse("L", no_case[named_strings]);
            BOOST_TEST(result);
            BOOST_TEST(*result == "50");
        }
        {
            auto const result = parse("l", no_case[named_strings]);
            BOOST_TEST(result);
            BOOST_TEST(*result == "50");
        }
    }
    // with mutation
    {
        symbols<int> roman_numerals;
        roman_numerals.insert_for_next_parse("I", 1);
        roman_numerals.insert_for_next_parse("V", 5);
        roman_numerals.insert_for_next_parse("X", 10);
        auto const add_numeral = [&roman_numerals](auto & context) {
            using namespace boost::parser::literals;
            char chars[2] = {get(_attr(context), 0_c), 0};
            roman_numerals.insert(context, chars, get(_attr(context), 1_c));
        };
        auto const numerals_parser = omit[roman_numerals] >>
                                     (char_ >> int_)[add_numeral] >>
                                     no_case[roman_numerals];

        {
            auto const result = parse("VL50L", numerals_parser);
            BOOST_TEST(result);
            BOOST_TEST(*result == 50);
            BOOST_TEST(!parse("L", roman_numerals));
        }
        {
            auto const result = parse("VL50l", numerals_parser);
            BOOST_TEST(result);
            BOOST_TEST(*result == 50);
            BOOST_TEST(!parse("L", roman_numerals));
        }
        {
            auto const result = parse("VC100C", numerals_parser);
            BOOST_TEST(result);
            BOOST_TEST(*result == 100);
            BOOST_TEST(!parse("C", roman_numerals));
        }
        {
            auto const result = parse("Vc100C", numerals_parser);
            BOOST_TEST(result);
            BOOST_TEST(*result == 100);
            BOOST_TEST(!parse("C", roman_numerals));
        }
    }
}

// multi_code_point_mapping)
{
    {
        constexpr auto capital_sharp_s_p = no_case[string(capital_sharp_s)];

        {
            auto const result =
                parse(null_term(capital_sharp_s) | as_utf32, capital_sharp_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)capital_sharp_s);
        }
        {
            auto const result =
                parse(null_term(small_sharp_s) | as_utf32, capital_sharp_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)small_sharp_s);
        }
        {
            auto const result = parse(null_term(double_s) | as_utf32, capital_sharp_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)double_s);
        }
        {
            BOOST_TEST(!parse("s" | as_utf32, capital_sharp_s_p));
        }
        {
            BOOST_TEST(!parse("sx" | as_utf32, capital_sharp_s_p));
        }
        {
            BOOST_TEST(!parse("xs" | as_utf32, capital_sharp_s_p));
        }
    }
    {
        constexpr auto small_sharp_s_p = no_case[string(small_sharp_s)];

        {
            auto const result =
                parse(null_term(capital_sharp_s) | as_utf32, small_sharp_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)capital_sharp_s);
        }
        {
            auto const result =
                parse(null_term(small_sharp_s) | as_utf32, small_sharp_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)small_sharp_s);
        }
        {
            auto const result =
                parse(null_term(double_s) | as_utf32, small_sharp_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)double_s);
        }
        {
            BOOST_TEST(!parse("s" | as_utf32, small_sharp_s_p));
        }
        {
            BOOST_TEST(!parse("sx" | as_utf32, small_sharp_s_p));
        }
        {
            BOOST_TEST(!parse("xs" | as_utf32, small_sharp_s_p));
        }
    }
    {
        constexpr auto double_s_p = no_case[string(double_s)];

        {
            auto const result =
                parse(null_term(capital_sharp_s) | as_utf32, double_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)capital_sharp_s);
        }
        {
            auto const result =
                parse(null_term(small_sharp_s) | as_utf32, double_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)small_sharp_s);
        }
        {
            auto const result =
                parse(null_term(double_s) | as_utf32, double_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)double_s);
        }
        {
            BOOST_TEST(!parse("s" | as_utf32, double_s_p));
        }
        {
            BOOST_TEST(!parse("sx" | as_utf32, double_s_p));
        }
        {
            BOOST_TEST(!parse("xs" | as_utf32, double_s_p));
        }
    }
    {
        constexpr auto s_p = no_case[string("s")];

        {
            auto const result = parse("s" | as_utf32, s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == "s");
        }
        {
            BOOST_TEST(!parse(null_term(capital_sharp_s) | as_utf32, s_p));
        }
        {
            BOOST_TEST(!parse(null_term(small_sharp_s) | as_utf32, s_p));
        }
        {
            BOOST_TEST(!parse(null_term(double_s) | as_utf32, s_p));
        }
        {
            BOOST_TEST(!parse("x" | as_utf32, s_p));
        }
    }
}

// detail_no_case_iter)
{
    {
        constexpr auto mixed_sharp_s1 = U"ẞs";
        std::string folded;
        auto first =
            detail::no_case_iter(mixed_sharp_s1, detail::text::null_sentinel);
        while (first != detail::text::null_sentinel) {
            folded.push_back((char)*first);
            ++first;
        }
        BOOST_TEST(folded == "sss");
        BOOST_TEST(first.base() == detail::text::null_sentinel);
    }
    {
        constexpr auto mixed_sharp_s2 = U"sẞ";
        std::string folded;
        auto first =
            detail::no_case_iter(mixed_sharp_s2, detail::text::null_sentinel);
        while (first != detail::text::null_sentinel) {
            folded.push_back((char)*first);
            ++first;
        }
        BOOST_TEST(folded == "sss");
        BOOST_TEST(first.base() == detail::text::null_sentinel);
    }
    {
        auto const street = U"Straße";
        std::string folded;
        auto const first_const =
            detail::no_case_iter(street, detail::text::null_sentinel);
        auto first = first_const;
        while (first != detail::text::null_sentinel) {
            folded.push_back((char)*first);
            ++first;
        }
        BOOST_TEST(folded == "strasse");
        BOOST_TEST(first.base() == detail::text::null_sentinel);

        first = first_const;
        std::u32string_view const sv = U"strasse";
        auto mismatches = detail::text::mismatch(
            first, detail::text::null_sentinel, sv.begin(), sv.end());
        BOOST_TEST(mismatches.first == detail::text::null_sentinel);
        BOOST_TEST(mismatches.second == sv.end());

        {
            first = first_const;
            auto search_result = detail::text::search(
                first, detail::text::null_sentinel, sv.begin(), sv.end());
            BOOST_TEST(search_result.begin() == first);
            BOOST_TEST(search_result.end() == detail::text::null_sentinel);
        }

        {
            first = first_const;
            auto search_result = detail::text::search(
                sv.begin(), sv.end(), first, detail::text::null_sentinel);
            BOOST_TEST(search_result.begin() == sv.begin());
            BOOST_TEST(search_result.end() == sv.end());
        }

        {
            detail::case_fold_array_t folded_char;
            auto folded_last = detail::case_fold('X', folded_char.begin());
            auto search_result = detail::text::search(
                sv.begin(), sv.end(), folded_char.begin(), folded_last);
            BOOST_TEST(search_result.begin() == sv.end());
            BOOST_TEST(search_result.end() == sv.end());
        }
    }
}

// detail_no_case_mismatch)
{
    constexpr auto mixed_sharp_s1 = U"ẞs";
    constexpr auto mixed_sharp_s2 = U"sẞ";
    auto const result = detail::no_case_aware_string_mismatch(
        mixed_sharp_s1,
        detail::text::null_sentinel,
        mixed_sharp_s2,
        detail::text::null_sentinel,
        true);
    BOOST_TEST(result.first == detail::text::null_sentinel);
    BOOST_TEST(result.second == detail::text::null_sentinel);
}

// longer_multi_code_point_mapping)
{
    constexpr auto mixed_sharp_s1 = u8"ẞs";
    constexpr auto mixed_sharp_s2 = u8"sẞ";

    constexpr auto mixed_sharp_s3 = u8"sßs";

    constexpr auto all_sharp_s1 = u8"ẞßss";
    constexpr auto all_sharp_s2 = u8"ssẞß";
    constexpr auto all_sharp_s3 = u8"ẞssß";

    constexpr auto triple_s = u8"sss";
    constexpr auto quadruple_s = u8"ssSs";

    {
        constexpr auto mixed_sharp_s1_p = no_case[string(mixed_sharp_s1)];
        {
            auto const result =
                parse(null_term(triple_s) | as_utf32, mixed_sharp_s1_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)triple_s);
        }
        {
            auto const result =
                parse(null_term(mixed_sharp_s1) | as_utf32, mixed_sharp_s1_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s1);
        }
        {
            auto const result =
                parse(null_term(mixed_sharp_s2) | as_utf32, mixed_sharp_s1_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s2);
        }
    }
    {
        constexpr auto mixed_sharp_s2_p = no_case[string(mixed_sharp_s2)];
        {
            auto const result =
                parse(null_term(triple_s) | as_utf32, mixed_sharp_s2_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)triple_s);
        }
        {
            auto const result =
                parse(null_term(mixed_sharp_s1) | as_utf32, mixed_sharp_s2_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s1);
        }
        {
            auto const result =
                parse(null_term(mixed_sharp_s2) | as_utf32, mixed_sharp_s2_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s2);
        }
    }
    {
        constexpr auto triple_s_p = no_case[string(triple_s)];
        {
            auto const result =
                parse(null_term(triple_s) | as_utf32, triple_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)triple_s);
        }
        {
            auto const result =
                parse(null_term(mixed_sharp_s1) | as_utf32, triple_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s1);
        }
        {
            auto const result =
                parse(null_term(mixed_sharp_s2) | as_utf32, triple_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s2);
        }
    }
    {
        constexpr auto mixed_sharp_s3_p = no_case[string(mixed_sharp_s3)];
        {
            auto const result =
                parse(null_term(mixed_sharp_s3) | as_utf32, mixed_sharp_s3_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s3);
        }
        {
            auto const result =
                parse(null_term(quadruple_s) | as_utf32, mixed_sharp_s3_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)quadruple_s);
        }
    }
    {
        constexpr auto quadruple_s_p = no_case[string(quadruple_s)];
        {
            auto const result =
                parse(null_term(mixed_sharp_s3) | as_utf32, quadruple_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)mixed_sharp_s3);
        }
        {
            auto const result =
                parse(null_term(quadruple_s) | as_utf32, quadruple_s_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)quadruple_s);
        }
    }
    {
        constexpr auto all_sharp_s1_p = no_case[string(all_sharp_s1)];
        {
            auto const result =
                parse(null_term(all_sharp_s1) | as_utf32, all_sharp_s1_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s1);
        }
        {
            auto const result =
                parse(null_term(all_sharp_s2) | as_utf32, all_sharp_s1_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s2);
        }
        {
            auto const result =
                parse(null_term(all_sharp_s3) | as_utf32, all_sharp_s1_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s3);
        }
        {
            auto const result = parse("ssssss" | as_utf32, all_sharp_s1_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == "ssssss");
        }
    }
    {
        constexpr auto all_sharp_s2_p = no_case[string(all_sharp_s1)];
        {
            auto const result =
                parse(null_term(all_sharp_s1) | as_utf32, all_sharp_s2_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s1);
        }
        {
            auto const result =
                parse(null_term(all_sharp_s2) | as_utf32, all_sharp_s2_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s2);
        }
        {
            auto const result =
                parse(null_term(all_sharp_s3) | as_utf32, all_sharp_s2_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s3);
        }
        {
            auto const result = parse("ssssss" | as_utf32, all_sharp_s2_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == "ssssss");
        }
    }
    {
        constexpr auto all_sharp_s3_p = no_case[string(all_sharp_s1)];
        {
            auto const result =
                parse(null_term(all_sharp_s1) | as_utf32, all_sharp_s3_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s1);
        }
        {
            auto const result =
                parse(null_term(all_sharp_s2) | as_utf32, all_sharp_s3_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s2);
        }
        {
            auto const result =
                parse(null_term(all_sharp_s3) | as_utf32, all_sharp_s3_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == (char const *)all_sharp_s3);
        }
        {
            auto const result = parse("ssSsss" | as_utf32, all_sharp_s3_p);
            BOOST_TEST(result);
            BOOST_TEST(*result == "ssSsss");
        }
    }
}

return boost::report_errors();
}
