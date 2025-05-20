/**
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/transform_replace.hpp>

#include <boost/core/lightweight_test.hpp>

#include "ill_formed.hpp"

#include <list>

#if (!defined(_MSC_VER) || BOOST_PARSER_USE_CONCEPTS) &&                       \
    (!defined(__GNUC__) || 12 <= __GNUC__ || !BOOST_PARSER_USE_CONCEPTS)

namespace bp = boost::parser;

auto f_str = [](std::vector<int> const & ints) {
    std::string retval;
    for (auto x : ints) {
        retval += std::to_string(x);
        retval += '_';
    }
    return retval;
};
auto f_str_ref = [](std::vector<int> const & ints) -> std::string & {
    static std::string retval;
    for (auto x : ints) {
        retval += std::to_string(x);
        retval += '_';
    }
    return retval;
};

#if BOOST_PARSER_USE_CONCEPTS
namespace deduction {
    using namespace std::literals;
    std::string str;
    auto const parser = bp::int_ % ',';
    auto const skip = bp::ws;
    using attr_t = std::vector<int>;

    auto deduced_1 = bp::transform_replace_view(
        str,
        parser,
        skip,
        bp::detail::utf_rvalue_shim<std::string, decltype(f_str), attr_t>(
            f_str),
        bp::trace::on);
    auto deduced_2 = bp::transform_replace_view(
        str,
        parser,
        skip,
        bp::detail::utf_rvalue_shim<std::string, decltype(f_str), attr_t>(
            f_str));
    auto deduced_3 = bp::transform_replace_view(
        str,
        parser,
        bp::detail::utf_rvalue_shim<std::string, decltype(f_str), attr_t>(
            f_str),
        bp::trace::on);
    auto deduced_4 = bp::transform_replace_view(
        str,
        parser,
        bp::detail::utf_rvalue_shim<std::string, decltype(f_str), attr_t>(
            f_str));
}
#endif

namespace detail_attr_type {
    constexpr auto int_char_p = bp::int_ >> bp::char_;

    static_assert(
        std::is_same_v<
            bp::detail::range_attr_t<std::string, decltype(int_char_p.parser_)>,
            bp::tuple<int, char>>);

    static_assert(
        std::is_same_v<
            bp::detail::
                range_attr_t<std::u32string, decltype(int_char_p.parser_)>,
            bp::tuple<int, char32_t>>);

    constexpr auto ints_p = *bp::int_;
    static_assert(
        std::is_same_v<
            bp::detail::range_attr_t<std::u32string, decltype(ints_p.parser_)>,
            std::vector<int>>);
}

#if defined(__cpp_char8_t)
auto f_u8str = [](std::vector<int> ints) {
    std::u8string retval;
    for (auto x : ints) {
        auto const s = std::to_string(x);
        retval.insert(retval.end(), s.begin(), s.end());
        retval += '_';
    }
    return retval;
};
// NOTE: *const* & return type!
auto f_u8str_ref = [](std::vector<int> ints) -> std::u8string const & {
    static std::u8string retval;
    for (auto x : ints) {
        auto const s = std::to_string(x);
        retval.insert(retval.end(), s.begin(), s.end());
        retval += '_';
    }
    return retval;
};
#endif
auto f_u16str = [](std::vector<int> ints) {
    std::u16string retval;
    for (auto x : ints) {
        auto const s = std::to_string(x);
        retval.insert(retval.end(), s.begin(), s.end());
        retval += '_';
    }
    return retval;
};
auto f_u16str_ref = [](std::vector<int> ints) -> std::u16string & {
    static std::u16string retval;
    for (auto x : ints) {
        auto const s = std::to_string(x);
        retval.insert(retval.end(), s.begin(), s.end());
        retval += '_';
    }
    return retval;
};
auto f_u32str = [](std::vector<int> ints) {
    std::u32string retval;
    for (auto x : ints) {
        auto const s = std::to_string(x);
        retval.insert(retval.end(), s.begin(), s.end());
        retval += '_';
    }
    return retval;
};
auto f_u32str_ref = [](std::vector<int> ints) -> std::u32string & {
    static std::u32string retval;
    for (auto x : ints) {
        auto const s = std::to_string(x);
        retval.insert(retval.end(), s.begin(), s.end());
        retval += '_';
    }
    return retval;
};

namespace detail_utf_rvalue_shim {
    constexpr auto ints_p = *bp::int_;

    using attr_t = std::vector<int>;

    // char -> char

    bp::detail::utf_rvalue_shim<std::string, decltype(f_str), attr_t>
        char_char_shim(f_str);
    static_assert(
        std::is_same_v<decltype(char_char_shim(attr_t{})), std::string &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(char_char_shim),
                  std::string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::string, decltype(f_str_ref), attr_t>
        char_char_ref_shim(f_str_ref);
    static_assert(
        std::is_same_v<decltype(char_char_ref_shim(attr_t{})), std::string &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(char_char_ref_shim),
                  std::string,
                  decltype(ints_p.parser_)>);

#if defined(__cpp_char8_t) && BOOST_PARSER_USE_CONCEPTS
    // char8_t -> char8_t

    bp::detail::utf_rvalue_shim<std::u8string, decltype(f_u8str), attr_t>
        u8_u8_shim(f_u8str);
    static_assert(
        std::is_same_v<decltype(u8_u8_shim(attr_t{})), std::u8string &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u8_u8_shim),
                  std::u8string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u8string, decltype(f_u8str_ref), attr_t>
        u8_u8_ref_shim(f_u8str_ref);
    static_assert(std::is_same_v<
                  decltype(u8_u8_ref_shim(attr_t{})),
                  std::u8string const &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u8_u8_ref_shim),
                  std::u8string,
                  decltype(ints_p.parser_)>);

    // char8_t -> char16_t

    bp::detail::utf_rvalue_shim<std::u8string, decltype(f_u16str), attr_t>
        u8_u16_shim(f_u16str);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u8_u16_shim),
                  std::u8string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u8string, decltype(f_u16str_ref), attr_t>
        u8_u16_ref_shim(f_u16str_ref);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u8_u16_ref_shim),
                  std::u8string,
                  decltype(ints_p.parser_)>);

    // char8_t -> char32_t

    bp::detail::utf_rvalue_shim<std::u8string, decltype(f_u32str), attr_t>
        u8_u32_shim(f_u32str);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u8_u32_shim),
                  std::u8string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u8string, decltype(f_u32str_ref), attr_t>
        u8_u32_ref_shim(f_u32str_ref);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u8_u32_ref_shim),
                  std::u8string,
                  decltype(ints_p.parser_)>);

    // char16_t -> char8_t

    bp::detail::utf_rvalue_shim<std::u16string, decltype(f_u8str), attr_t>
        u16_u8_shim(f_u8str);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u16_u8_shim),
                  std::u16string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u16string, decltype(f_u8str_ref), attr_t>
        u16_u8_ref_shim(f_u8str_ref);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u16_u8_ref_shim),
                  std::u16string,
                  decltype(ints_p.parser_)>);

    // char32_t -> char8_t

    bp::detail::utf_rvalue_shim<std::u32string, decltype(f_u8str), attr_t>
        u32_u8_shim(f_u8str);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u32_u8_shim),
                  std::u32string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u32string, decltype(f_u8str_ref), attr_t>
        u32_u8_ref_shim(f_u8str_ref);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u32_u8_ref_shim),
                  std::u32string,
                  decltype(ints_p.parser_)>);
#endif
    // char16_t -> char16_t

    bp::detail::utf_rvalue_shim<std::u16string, decltype(f_u16str), attr_t>
        u16_u16_shim(f_u16str);
    static_assert(
        std::is_same_v<decltype(u16_u16_shim(attr_t{})), std::u16string &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u16_u16_shim),
                  std::u16string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u16string, decltype(f_u16str_ref), attr_t>
        u16_u16_ref_shim(f_u16str_ref);
    static_assert(
        std::is_same_v<decltype(u16_u16_ref_shim(attr_t{})), std::u16string &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u16_u16_ref_shim),
                  std::u16string,
                  decltype(ints_p.parser_)>);

    // char16_t -> char32_t

    bp::detail::utf_rvalue_shim<std::u16string, decltype(f_u32str), attr_t>
    u16_u32_shim(f_u32str);
#if !BOOST_PARSER_USE_CONCEPTS
    static_assert(
        std::is_same_v<
            decltype(u16_u32_shim(attr_t{})),
            bp::detail::text::utf16_view<
                bp::detail::text::detail::owning_view<std::u32string>> &>);
#endif
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u16_u32_shim),
                  std::u16string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u16string, decltype(f_u32str_ref), attr_t>
        u16_u32_ref_shim(f_u32str_ref);
#if !BOOST_PARSER_USE_CONCEPTS
    static_assert(std::is_same_v<
                  decltype(u16_u32_ref_shim(attr_t{})),
                  bp::detail::text::utf16_view<
                      bp::detail::text::detail::ref_view<std::u32string>> &>);
#endif
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u16_u32_ref_shim),
                  std::u16string,
                  decltype(ints_p.parser_)>);

    // char32_t -> char32_t

    bp::detail::utf_rvalue_shim<std::u32string, decltype(f_u32str), attr_t>
        u32_u32_shim(f_u32str);
    static_assert(
        std::is_same_v<decltype(u32_u32_shim(attr_t{})), std::u32string &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u32_u32_shim),
                  std::u32string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u32string, decltype(f_u32str_ref), attr_t>
        u32_u32_ref_shim(f_u32str_ref);
    static_assert(
        std::is_same_v<decltype(u32_u32_ref_shim(attr_t{})), std::u32string &>);
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u32_u32_ref_shim),
                  std::u32string,
                  decltype(ints_p.parser_)>);

    // char32_t -> char16_t

    bp::detail::utf_rvalue_shim<std::u32string, decltype(f_u16str), attr_t>
        u32_u16_shim(f_u16str);
#if !BOOST_PARSER_USE_CONCEPTS
    static_assert(
        std::is_same_v<
            decltype(u32_u16_shim(attr_t{})),
            bp::detail::text::utf32_view<
                bp::detail::text::detail::owning_view<std::u16string>> &>);
#endif
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u32_u16_shim),
                  std::u32string,
                  decltype(ints_p.parser_)>);

    bp::detail::utf_rvalue_shim<std::u32string, decltype(f_u16str_ref), attr_t>
        u32_u16_ref_shim(f_u16str_ref);
#if !BOOST_PARSER_USE_CONCEPTS
    static_assert(std::is_same_v<
                  decltype(u32_u16_ref_shim(attr_t{})),
                  bp::detail::text::utf32_view<
                      bp::detail::text::detail::ref_view<std::u16string>> &>);
#endif
    static_assert(bp::detail::transform_replacement_for<
                  decltype(u32_u16_ref_shim),
                  std::u32string,
                  decltype(ints_p.parser_)>);
}

int main()
{

// detail_attr_search_repack_shim
{
    using namespace bp::literals;

    {
        std::string str = "";
        auto parser = bp::string("XYZ");

        // Follows body of attr_search_impl() that constructs a custom parser
        // from the given one.
        auto first = bp::detail::text::detail::begin(str);
        auto const last = bp::detail::text::detail::end(str);
        auto match_first = first;
        auto match_last = first;
        auto before = [&match_first](auto & ctx) {
            match_first = _where(ctx).begin();
        };
        auto after = [&match_last](auto & ctx) {
            match_last = _where(ctx).begin();
        };
        auto const search_parser =
            bp::omit[*(bp::char_ - parser)] >>
            -bp::lexeme[bp::eps[before] >> bp::skip[parser] >> bp::eps[after]];

        auto result = bp::prefix_parse(
            first, last, search_parser, bp::ws, bp::trace::off);
        static_assert(std::is_same_v<
                      decltype(result),
                      std::optional<std::optional<std::string>>>);
        static_assert(std::is_same_v<
                      decltype(bp::prefix_parse(
                          first, last, search_parser, bp::ws, bp::trace::off)),
                      std::optional<std::optional<std::string>>>);
    }
    {
        std::string str = "";
        auto parser = bp::string("XYZ");
        bp::detail::attr_search_impl(str, parser, bp::ws, bp::trace::off);
    }
    {
        std::string str = "";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == std::begin(str));
        BOOST_TEST(subrng.end() == std::begin(str));
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "");
    }
    {
        char const str[] = "not here";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == std::end(str));
        BOOST_TEST(subrng.end() == std::end(str));
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "");
    }
    {
        char const str[] = "aaXYZb";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == str + 2);
        BOOST_TEST(subrng.end() == str + 5);
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "XYZ");
    }
    {
        char const str[] = "XYZab";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == str + 0);
        BOOST_TEST(subrng.end() == str + 3);
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "XYZ");
    }
    {
        char const str[] = "gbXYZ";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == str + 2);
        BOOST_TEST(subrng.end() == str + 5);
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "XYZ");
    }
    {
        char const str[] = "XYZ";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == str + 0);
        BOOST_TEST(subrng.end() == str + 3);
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "XYZ");
    }
    {
        char const str[] = "XXYZZ";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == str + 1);
        BOOST_TEST(subrng.end() == str + 4);
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "XYZ");
    }
    {
        char const str[] = "XXYZZ";
        auto result = bp::detail::attr_search_repack_shim(
            str, bp::string("XYZ"), bp::ws, bp::trace::off);
        auto subrng = bp::get(result, 0_c);
        BOOST_TEST(subrng.begin() == str + 1);
        BOOST_TEST(subrng.end() == str + 4);
        auto result_str = bp::get(result, 1_c);
        BOOST_TEST(result_str == "XYZ");
    }
}

// transform_replace
{
    {
        auto r = bp::transform_replace("", bp::int_ % ',', bp::ws, f_str);
        int count = 0;
        for (auto subrange : r) {
            (void)subrange;
            ++count;
        }
        BOOST_TEST(count == 0);
    }
    {
        char const str[] = "ab c 1,   2, 3  d e f";
        auto r = bp::transform_replace(str, bp::int_ % ',', bp::ws, f_str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(count == 3);
        BOOST_TEST(replace_result == "ab c 1_2_3_  d e f");
    }
    {
        char const str[] = "ab c 1,   2, 3  d e f";
        auto r = bp::transform_replace(str, bp::int_ % ',', bp::ws, f_str_ref);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(count == 3);
        BOOST_TEST(replace_result == "ab c 1_2_3_  d e f");
    }
    {
        char const str[] = "a a 1,2,3baa ba1 ,2 , 3";
        auto r = str | bp::transform_replace(
                           bp::int_ % ',', bp::ws, f_str, bp::trace::off);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "a a 1_2_3_baa ba1_2_3_");
        BOOST_TEST(count == 4);
    }
    {
        char const str[] = "aa1,2,3baaba1,2,3 4,5,6";
        auto r = str | bp::transform_replace(bp::int_ % ',', f_str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "aa1_2_3_baaba1_2_3_ 4_5_6_");
        BOOST_TEST(count == 6);
    }
    {
        char const str[] = "0,0aa1,2,3baaba1,2,3 4,5,6";
        auto r = str | bp::transform_replace(bp::int_ % ',', f_str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "0_0_aa1_2_3_baaba1_2_3_ 4_5_6_");
        BOOST_TEST(count == 7);
    }
    {
        char const str[] = "88,88 0,0aa1,2,3baaba1,2,3 4,5,6";
        auto r = str | bp::transform_replace(bp::int_ % ',', f_str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "88_88_ 0_0_aa1_2_3_baaba1_2_3_ 4_5_6_");
        BOOST_TEST(count == 9);
    }
}

// transform_replace_unicode
{
    {
        char const str_[] = "";
        auto str = str_ | bp::as_utf8;
        auto r = bp::transform_replace(str, bp::int_ % ',', bp::ws, f_u16str);
        int count = 0;
        for (auto subrange : r) {
            (void)subrange;
            ++count;
        }
        BOOST_TEST(count == 0);
    }
    {
        char const * str_ = "aa2,3,4b";
        auto str = bp::null_term(str_) | bp::as_utf16;
        auto r = bp::transform_replace(str, bp::int_ % ',', bp::ws, f_u16str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "aa2_3_4_b");
        BOOST_TEST(count == 3);
    }
    {
        char const str_[] = "a a 3,4,5 baaba7, 8 ,9";
        auto str = str_ | bp::as_utf32;
        auto r = str | bp::transform_replace(
                           bp::int_ % ',', bp::ws, f_u32str, bp::trace::off);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "a a 3_4_5_ baaba7_8_9_");
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "aa88,99baaba111,2222";
        auto str = str_ | bp::as_utf8;
        const auto r = str | bp::transform_replace(
                                 bp::int_ % ',', f_u16str, bp::trace::off);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "aa88_99_baaba111_2222_");
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "aa88,99baaba111,2222";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::transform_replace(bp::int_ % ',', f_u32str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "aa88_99_baaba111_2222_");
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "aa88,99baaba111,2222 3,4";
        auto str = str_ | bp::as_utf32;
        auto r = str | bp::transform_replace(bp::int_ % ',', f_u16str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "aa88_99_baaba111_2222_ 3_4_");
        BOOST_TEST(count == 6);
    }
    {
        char const str_[] = "1aa88,99baaba111,2222 3,4";
        auto str = str_ | bp::as_utf8;
        auto r = str | bp::transform_replace(bp::int_ % ',', f_u32str);
        int count = 0;
        std::string replace_result;
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            replace_result += str;
            ++count;
        }
        BOOST_TEST(replace_result == "1_aa88_99_baaba111_2222_ 3_4_");
        BOOST_TEST(count == 7);
    }
}

#if BOOST_PARSER_USE_CONCEPTS && (!defined(__GNUC__) || 12 <= __GNUC__)
// Older GCCs don't like the use of temporaries like the std::string("foo")
// below.  This causes | join to break.
// join_compat)
{
    {
        char const str_[] = "1aa88,99baaba111,2222 3,4";
        auto str = str_ | bp::as_utf16;
        auto rng = str | bp::transform_replace(bp::int_ % ',', f_u32str) |
                   std::views::join;
        std::string transform_replace_result;
        for (auto ch : rng) {
            static_assert(std::is_same_v<decltype(ch), char16_t>);
            transform_replace_result.push_back((char)ch);
        }
        BOOST_TEST(transform_replace_result == "1_aa88_99_baaba111_2222_ 3_4_");
    }
    {
        char const str[] = "1aa88,99baaba111,2222 3,4";
        auto rng = str | bp::as_utf32 |
                   bp::transform_replace(bp::int_ % ',', f_u8str) |
                   std::views::join;
        std::string transform_replace_result;
        for (auto ch : rng) {
            static_assert(std::is_same_v<decltype(ch), char32_t>);
            transform_replace_result.push_back((char)ch);
        }
        BOOST_TEST(transform_replace_result == "1_aa88_99_baaba111_2222_ 3_4_");
    }

    {
        char const str[] = "1aa88,99baaba111,2222 3,4";
        auto rng = str | bp::transform_replace(bp::int_ % ',', f_str) |
                   std::views::join;
        std::string transform_replace_result;
        for (auto ch : rng) {
            transform_replace_result.push_back(ch);
        }
        BOOST_TEST(transform_replace_result == "1_aa88_99_baaba111_2222_ 3_4_");
    }
    {
        std::string str = "1aa88,99baaba111,2222 3,4";
        auto rng = str | bp::transform_replace(bp::int_ % ',', f_str) |
                   std::views::join;
        std::string transform_replace_result;
        for (auto ch : rng) {
            transform_replace_result.push_back(ch);
        }
        BOOST_TEST(transform_replace_result == "1_aa88_99_baaba111_2222_ 3_4_");
    }
    {
        std::string const str = "1aa88,99baaba111,2222 3,4";
        auto rng = str | bp::transform_replace(bp::int_ % ',', f_str) |
                   std::views::join;
        std::string transform_replace_result;
        for (auto ch : rng) {
            transform_replace_result.push_back(ch);
        }
        BOOST_TEST(transform_replace_result == "1_aa88_99_baaba111_2222_ 3_4_");
    }
    {
        auto rng = std::string("1aa88,99baaba111,2222 3,4") |
                   bp::transform_replace(bp::int_ % ',', f_str) |
                   std::views::join;
        std::string transform_replace_result;
        for (auto ch : rng) {
            transform_replace_result.push_back(ch);
        }
        BOOST_TEST(transform_replace_result == "1_aa88_99_baaba111_2222_ 3_4_");
    }
}
#endif

// doc_examples
{
    {
        auto string_sum = [](std::vector<int> const & ints) {
            return std::to_string(std::accumulate(ints.begin(), ints.end(), 0));
        };

        auto rng = "There are groups of [1, 2, 3, 4, 5] in the set." |
                   bp::transform_replace(
                       '[' >> bp::int_ % ',' >> ']', bp::ws, string_sum);
        int count = 0;
        // Prints "There are groups of 15 in the set".
        for (auto subrange : rng) {
            for (auto ch : subrange) {
                std::cout << ch;
            }
            ++count;
        }
        std::cout << "\n";
        assert(count == 3);
    }
}

return boost::report_errors();
}

#else
int main() { return boost::report_errors(); }
#endif
