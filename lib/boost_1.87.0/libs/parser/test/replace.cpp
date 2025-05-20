/**
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/replace.hpp>

#include <boost/core/lightweight_test.hpp>

#include "ill_formed.hpp"

#include <list>

#if defined(_MSC_VER) && !BOOST_PARSER_USE_CONCEPTS
int main() {}
#else

namespace bp = boost::parser;

#if BOOST_PARSER_USE_CONCEPTS
namespace deduction {
    using namespace std::literals;
    std::string str;
    auto const parser = bp::char_;
    auto const skip = bp::ws;

    auto deduced_1 = bp::replace_view(str, parser, skip, "foo", bp::trace::on);
    auto deduced_2 = bp::replace_view(str, parser, skip, "foo");
    auto deduced_3 = bp::replace_view(str, parser, "foo", bp::trace::on);
    auto deduced_4 = bp::replace_view(str, parser, "foo");
}
#endif

// MSVC and older Clangs produce hard errors here, so ill_formed does not
// work.
#if defined(__cpp_char8_t) && !defined(_MSC_VER) &&                            \
    (!defined(__clang_major__) || 16 <= __clang_major__)
char const empty_str[] = "";

template<typename T>
using char_str_utf8_replacement =
    decltype(std::declval<T>() | bp::replace(bp::lit("XYZ"), std::declval<T>() | bp::as_utf8));
static_assert(ill_formed<char_str_utf8_replacement, decltype(empty_str)>{});

template<typename T>
using char_str_utf16_replacement =
    decltype(std::declval<T>() | bp::replace(bp::lit("XYZ"), std::declval<T>() | bp::as_utf16));
static_assert(ill_formed<char_str_utf16_replacement, decltype(empty_str)>{});
#endif

static_assert(
    bp::detail::range_utf_format<char const *&>() == bp::detail::no_format);

int main()
{

// bind_back
{
    auto f = [](auto x) {
        static_assert(std::is_same_v<
                      decltype(x),
                      BOOST_PARSER_SUBRANGE<char const *, char const *>>);
        BOOST_TEST(x.size() == 4u); // stripped off null
    };
    bp::detail::stl_interfaces::bind_back(f, "text");
}

// either_iterator
{
    {
        std::list<int> l({1, 2, 3});
        std::vector<int> v({4, 5, 6});
        bp::detail::either_iterator<std::list<int>, std::vector<int>>
            either_l_begin(l.begin());
        bp::detail::either_iterator<std::list<int>, std::vector<int>>
            either_l_end(l.end());
        bp::detail::either_iterator<std::list<int>, std::vector<int>>
            either_v_begin(v.begin());
        bp::detail::either_iterator<std::list<int>, std::vector<int>>
            either_v_end(v.end());

        int const l_array[] = {1, 2, 3};
        auto l_array_curr = l_array;
        for (auto it = either_l_begin; it != either_l_end;
             ++it, ++l_array_curr) {
            BOOST_TEST(*it == *l_array_curr);
        }

        int const v_array[] = {4, 5, 6};
        auto v_array_curr = v_array;
        for (auto it = either_v_begin; it != either_v_end;
             ++it, ++v_array_curr) {
            BOOST_TEST(*it == *v_array_curr);
        }
    }
}

// replace
{
    {
        auto r = bp::replace("", bp::lit("XYZ"), bp::ws, "foo");
        int count = 0;
        for (auto subrange : r) {
            (void)subrange;
            ++count;
        }
        BOOST_TEST(count == 0);
    }
    {
        char const str[] = "aaXYZb";
        auto r = bp::replace(str, bp::lit("XYZ"), bp::ws, "foo");
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "b"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str[] = "a a XYZ baa ba XYZ";
        auto r =
            str | bp::replace(bp::lit("XYZ"), bp::ws, "foo", bp::trace::off);
        int count = 0;
        std::string_view const strs[] = {"a a ", "foo", " baa ba ", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
#if !defined(__GNUC__) || 12 <= __GNUC__
    // Older GCCs don't like the use of temporaries like the
    // std::string("foo") below.
    {
        char const str[] = "aaXYZbaabaXYZ";
        auto r = str | bp::replace(
                           bp::lit("XYZ"), std::string("foo"), bp::trace::off);
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "baaba", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
#endif
    {
        char const str[] = "aaXYZbaabaXYZ";
        const auto r = str | bp::replace(bp::lit("XYZ"), "foo");
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "baaba", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str[] = "aaXYZbaabaXYZXYZ";
        auto r = str | bp::replace(bp::lit("XYZ"), "foo");
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
    {
        char const str[] = "XYZaaXYZbaabaXYZXYZ";
        auto r = str | bp::replace(bp::lit("XYZ"), "foo");
        int count = 0;
        std::string_view const strs[] = {
            "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 6);
    }
    {
        char const str[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto r = str | bp::replace(bp::lit("XYZ"), "foo");
        int count = 0;
        std::string_view const strs[] = {
            "foo", "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 7);
    }
    {
        char const * str = "XYZXYZaaXYZbaabaXYZXYZ";
        char const * replacement = "foo";
        auto r = bp::null_term(str) |
            bp::replace(bp::lit("XYZ"), bp::null_term(replacement));
        int count = 0;
        std::string_view const strs[] = {
            "foo", "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 7);
    }
    {
        char const * str = "XYZXYZaaXYZbaabaXYZXYZ";
        char const * replacement = "foo";
        auto const r = bp::null_term(str) |
            bp::replace(bp::lit("XYZ"), bp::null_term(replacement));
        int count = 0;
        std::string_view const strs[] = {
            "foo", "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 7);
    }
}

// replace_unicode
{
    {
        char const str_[] = "";
        auto str = str_ | bp::as_utf8;
        auto r = bp::replace(str, bp::lit("XYZ"), bp::ws, "foo" | bp::as_utf8);
        int count = 0;
        for (auto subrange : r) {
            (void)subrange;
            ++count;
        }
        BOOST_TEST(count == 0);
    }
    {
        char const * str_ = "aaXYZb";
        auto str = bp::null_term(str_) | bp::as_utf16;
        auto r = bp::replace(str, bp::lit("XYZ"), bp::ws, "foo" | bp::as_utf16);
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "b"};
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf32;
        auto r =
            str |
            bp::replace(
                bp::lit("XYZ"), bp::ws, "foo" | bp::as_utf32, bp::trace::off);
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "baaba", "foo"};
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf8;
        auto r = str | bp::replace(
                           bp::lit("XYZ"), "foo" | bp::as_utf8, bp::trace::off);
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "baaba", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::replace(bp::lit("XYZ"), "foo");
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "baaba", "foo"};
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "aaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf32;
        auto r = str | bp::replace(bp::lit("XYZ"), "foo");
        int count = 0;
        std::string_view const strs[] = {"aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
    {
        char const str_[] = "XYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf8;
        auto r = str | bp::replace(bp::lit("XYZ"), "foo" | bp::as_utf8);
        int count = 0;
        std::string_view const strs[] = {
            "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            std::string str(subrange.begin(), subrange.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 6);
    }
    {
        char const str_[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::replace(bp::lit("XYZ"), "foo");
        int count = 0;
        std::string_view const strs[] = {
            "foo", "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 7);
    }
    {
        char const str_[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::replace(bp::lit("XYZ"), "foo" | bp::as_utf8);
        int count = 0;
        std::string_view const strs[] = {
            "foo", "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 7);
    }
    {
        char const str_[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::replace(bp::lit("XYZ"), "foo" | bp::as_utf32);
        int count = 0;
        std::string_view const strs[] = {
            "foo", "foo", "aa", "foo", "baaba", "foo", "foo"};
        for (auto subrange : r) {
            auto u8sub = subrange | bp::as_utf8;
            std::string str(u8sub.begin(), u8sub.end());
            BOOST_TEST(str == strs[count]);
            ++count;
        }
        BOOST_TEST(count == 7);
    }
}

#if BOOST_PARSER_USE_CONCEPTS && (!defined(__GNUC__) || 12 <= __GNUC__)
// Older GCCs don't like the use of temporaries like the std::string("foo")
// below.  This causes | join to break.
// join_compat)
{
    {
        char const str[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto rng = str | bp::as_utf32 |
                   bp::replace(bp::lit("XYZ"), "foo" | bp::as_utf8) |
                   std::views::join;
        std::string replace_result;
        for (auto ch : rng) {
            static_assert(std::is_same_v<decltype(ch), char32_t>);
            replace_result.push_back((char)ch);
        }
        BOOST_TEST(replace_result == "foofooaafoobaabafoofoo");
    }

    {
        char const str[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto rng = str | bp::replace(bp::lit("XYZ"), "foo") | std::views::join;
        std::string replace_result;
        for (auto ch : rng) {
            replace_result.push_back(ch);
        }
        BOOST_TEST(replace_result == "foofooaafoobaabafoofoo");
    }
    {
        std::string str = "XYZXYZaaXYZbaabaXYZXYZ";
        auto rng = str | bp::replace(bp::lit("XYZ"), "foo") | std::views::join;
        std::string replace_result;
        for (auto ch : rng) {
            replace_result.push_back(ch);
        }
        BOOST_TEST(replace_result == "foofooaafoobaabafoofoo");
    }
    {
        std::string const str = "XYZXYZaaXYZbaabaXYZXYZ";
        auto rng = str | bp::replace(bp::lit("XYZ"), "foo") | std::views::join;
        std::string replace_result;
        for (auto ch : rng) {
            replace_result.push_back(ch);
        }
        BOOST_TEST(replace_result == "foofooaafoobaabafoofoo");
    }
    {
        auto rng = std::string("XYZXYZaaXYZbaabaXYZXYZ") |
                   bp::replace(bp::lit("XYZ"), "foo") | std::views::join;
        std::string replace_result;
        for (auto ch : rng) {
            replace_result.push_back(ch);
        }
        BOOST_TEST(replace_result == "foofooaafoobaabafoofoo");
    }
}
#endif

// doc_examples
{
    // clang-format off
    {
        namespace bp = boost::parser;
        auto card_number = bp::int_ >> bp::repeat(3)['-' >> bp::int_];
        auto rng = "My credit card number is 1234-5678-9012-3456." | bp::replace(card_number, "XXXX-XXXX-XXXX-XXXX");
        int count = 0;
        // Prints My credit card number is XXXX-XXXX-XXXX-XXXX.
        for (auto subrange : rng) {
            std::cout << std::string_view(subrange.begin(), subrange.end() - subrange.begin());
            ++count;
        }
        std::cout << "\n";
        assert(count == 3);
    }
#if BOOST_PARSER_USE_CONCEPTS && (!defined(__GNUC__) || 12 <= __GNUC__)
    {
        namespace bp = boost::parser;
        auto card_number = bp::int_ >> bp::repeat(3)['-' >> bp::int_];
        auto rng = "My credit card number is 1234-5678-9012-3456." |
                   bp::replace(card_number, "XXXX-XXXX-XXXX-XXXX") |
                   std::views::join;
        std::string replace_result;
        for (auto ch : rng) {
            replace_result.push_back(ch);
        }
        assert(replace_result == "My credit card number is XXXX-XXXX-XXXX-XXXX.");
    }
#endif
    // clang-format on
}

    return boost::report_errors();
}

#endif
