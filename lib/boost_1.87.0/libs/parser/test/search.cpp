/**
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/search.hpp>

#include <boost/core/lightweight_test.hpp>


namespace bp = boost::parser;

#if BOOST_PARSER_USE_CONCEPTS
namespace deduction {
    std::string str;
    auto const parser = bp::char_;
    auto const skip = bp::ws;

    auto deduced_1 = bp::search_all_view(str, parser, skip, bp::trace::on);
    auto deduced_2 = bp::search_all_view(str, parser, skip);
    auto deduced_3 = bp::search_all_view(str, parser, bp::trace::on);
    auto deduced_4 = bp::search_all_view(str, parser);
}
#endif

int main()
{

// search_range_skip
{
    // array of char
    {
        char const str[] = "";
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str);
        BOOST_TEST(result.end() == str);
    }
    {
        char const str[] = "not here";
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::end(str) - 1);
        BOOST_TEST(result.end() == std::end(str) - 1);
    }
    {
        char const str[] = "aaXYZb";
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const str[] = "XYZab";
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "gbXYZ";
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const str[] = "XYZ";
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "XXYZZ";
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }
    {
        char const * str = "XXYZZ";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }
#if BOOST_PARSER_USE_CONCEPTS
    {
        auto result = bp::search(std::string("XXYZZ"), bp::lit("XYZ"), bp::ws);
        static_assert(std::same_as<decltype(result), std::ranges::dangling>);
    }
#endif

    // array of char | as_utf32
    {
        char const str[] = "";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str);
        BOOST_TEST(result.end() == str);
    }
    {
        char const str[] = "not here";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::end(str) - 1);
        BOOST_TEST(result.end() == std::end(str) - 1);
    }
    {
        char const str[] = "aaXYZb";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const str[] = "XYZab";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "gbXYZ";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const str[] = "XYZ";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "XXYZZ";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }

    // pointer
    {
        char const * str = "";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str);
        BOOST_TEST(result.end() == str);
    }
    {
        char const * str = "not here";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == bp::null_sentinel_t{});
        BOOST_TEST(result.end() == bp::null_sentinel_t{});
    }
    {
        char const * str = "aaXYZb";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const * str = "XYZab";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const * str = "gbXYZ";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const * str = "XYZ";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const * str = "XXYZZ";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }

    // pointer
    {
        char const * str_ = "";
        auto str = bp::null_term(str_) | bp::as_utf32;
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str.begin());
        BOOST_TEST(result.end() == str.end());
    }
    {
        char const * str_ = "not here";
        auto str = bp::null_term(str_) | bp::as_utf16;
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str.end());
        BOOST_TEST(result.end() == str.end());
    }
    {
        char const * str_ = "aaXYZb";
        auto str = bp::null_term(str_) | bp::as_utf8;
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 2));
        BOOST_TEST(result.end() == std::next(str.begin(), 5));
    }
    {
        char const * str_ = "XYZab";
        auto str = bp::null_term(str_) | bp::as_utf32;
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
    {
        char const * str_ = "gbXYZ";
        auto str = bp::null_term(str_) | bp::as_utf16;
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 2));
        BOOST_TEST(result.end() == std::next(str.begin(), 5));
    }
    {
        char const * str_ = "XYZ";
        auto str = bp::null_term(str_) | bp::as_utf8;
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
    {
        char const * str_ = "XXYZZ";
        auto str = bp::null_term(str_) | bp::as_utf32;
        auto result = bp::search(str, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 1));
        BOOST_TEST(result.end() == std::next(str.begin(), 4));
    }
}

// search_iters_skip
{
    // array of char
    {
        char const str[] = "XYZab";
        auto result =
            bp::search(std::begin(str), std::end(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "gbXYZ";
        auto result =
            bp::search(std::begin(str), std::end(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const str[] = "XYZ";
        auto result =
            bp::search(std::begin(str), std::end(str), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }

    // array of char | as_utf32
    {
        char const str_[] = "";
        auto str = str_ | bp::as_utf32;
        auto result =
            bp::search(str.begin(), str.end(), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str_);
        BOOST_TEST(result.end() == str_);
    }
    {
        char const str_[] = "XYZ";
        auto str = str_ | bp::as_utf32;
        auto result =
            bp::search(str.begin(), str.end(), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str_ + 0);
        BOOST_TEST(result.end() == str_ + 3);
    }
    {
        char const str_[] = "XXYZZ";
        auto str = str_ | bp::as_utf32;
        auto result =
            bp::search(str.begin(), str.end(), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str_ + 1);
        BOOST_TEST(result.end() == str_ + 4);
    }

    // pointer
    {
        char const * str = "";
        auto result =
            bp::search(str, bp::null_sentinel_t{}, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str);
        BOOST_TEST(result.end() == str);
    }
    {
        char const * str = "not here";
        auto result =
            bp::search(str, bp::null_sentinel_t{}, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == bp::null_sentinel_t{});
        BOOST_TEST(result.end() == bp::null_sentinel_t{});
    }
    {
        char const * str = "XXYZZ";
        auto result =
            bp::search(str, bp::null_sentinel_t{}, bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }

    // pointer
    {
        char const * str_ = "XYZab";
        auto str = bp::null_term(str_) | bp::as_utf32;
        auto result =
            bp::search(str.begin(), str.end(), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
    {
        char const * str_ = "gbXYZ";
        auto str = bp::null_term(str_) | bp::as_utf16;
        auto result =
            bp::search(str.begin(), str.end(), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 2));
        BOOST_TEST(result.end() == std::next(str.begin(), 5));
    }
    {
        char const * str_ = "XYZ";
        auto str = bp::null_term(str_) | bp::as_utf8;
        auto result =
            bp::search(str.begin(), str.end(), bp::lit("XYZ"), bp::ws);
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
}

// search_range_no_skip
{
    // array of char
    {
        char const str[] = "XYZab";
        auto result = bp::search(str, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "gbXYZ";
        auto result = bp::search(str, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const str[] = "XYZ";
        auto result = bp::search(str, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }

    // array of char | as_utf32
    {
        char const str[] = "";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str);
        BOOST_TEST(result.end() == str);
    }
    {
        char const str[] = "XYZ";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "XXYZZ";
        auto result = bp::search(str | bp::as_utf32, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }

    // pointer
    {
        char const * str = "";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str);
        BOOST_TEST(result.end() == str);
    }
    {
        char const * str = "not here";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == bp::null_sentinel_t{});
        BOOST_TEST(result.end() == bp::null_sentinel_t{});
    }
    {
        char const * str = "XXYZZ";
        auto result = bp::search(bp::null_term(str), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }

    // pointer
    {
        char const * str_ = "XYZab";
        auto str = bp::null_term(str_) | bp::as_utf32;
        auto result = bp::search(str, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
    {
        char const * str_ = "gbXYZ";
        auto str = bp::null_term(str_) | bp::as_utf16;
        auto result = bp::search(str, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == std::next(str.begin(), 2));
        BOOST_TEST(result.end() == std::next(str.begin(), 5));
    }
    {
        char const * str_ = "XYZ";
        auto str = bp::null_term(str_) | bp::as_utf8;
        auto result = bp::search(str, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
}

// search_iters_no_skip
{
    // array of char
    {
        char const str[] = "XYZab";
        auto result =
            bp::search(std::begin(str), std::end(str), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }
    {
        char const str[] = "gbXYZ";
        auto result =
            bp::search(std::begin(str), std::end(str), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 2);
        BOOST_TEST(result.end() == str + 5);
    }
    {
        char const str[] = "XYZ";
        auto result =
            bp::search(std::begin(str), std::end(str), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 0);
        BOOST_TEST(result.end() == str + 3);
    }

    // array of char | as_utf32
    {
        char const str_[] = "";
        auto str = str_ | bp::as_utf32;
        auto result = bp::search(str.begin(), str.end(), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str_);
        BOOST_TEST(result.end() == str_);
    }
    {
        char const str_[] = "XYZ";
        auto str = str_ | bp::as_utf32;
        auto result = bp::search(str.begin(), str.end(), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str_ + 0);
        BOOST_TEST(result.end() == str_ + 3);
    }
    {
        char const str_[] = "XXYZZ";
        auto str = str_ | bp::as_utf32;
        auto result = bp::search(str.begin(), str.end(), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str_ + 1);
        BOOST_TEST(result.end() == str_ + 4);
    }

    // pointer
    {
        char const * str = "";
        auto result = bp::search(str, bp::null_sentinel_t{}, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str);
        BOOST_TEST(result.end() == str);
    }
    {
        char const * str = "not here";
        auto result = bp::search(str, bp::null_sentinel_t{}, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == bp::null_sentinel_t{});
        BOOST_TEST(result.end() == bp::null_sentinel_t{});
    }
    {
        char const * str = "XXYZZ";
        auto result = bp::search(str, bp::null_sentinel_t{}, bp::lit("XYZ"));
        BOOST_TEST(result.begin() == str + 1);
        BOOST_TEST(result.end() == str + 4);
    }

    // pointer
    {
        char const * str_ = "XYZab";
        auto str = bp::null_term(str_) | bp::as_utf32;
        auto result = bp::search(str.begin(), str.end(), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
    {
        char const * str_ = "gbXYZ";
        auto str = bp::null_term(str_) | bp::as_utf16;
        auto result = bp::search(str.begin(), str.end(), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == std::next(str.begin(), 2));
        BOOST_TEST(result.end() == std::next(str.begin(), 5));
    }
    {
        char const * str_ = "XYZ";
        auto str = bp::null_term(str_) | bp::as_utf8;
        auto result = bp::search(str.begin(), str.end(), bp::lit("XYZ"));
        BOOST_TEST(result.begin() == std::next(str.begin(), 0));
        BOOST_TEST(result.end() == std::next(str.begin(), 3));
    }
}

// search_all
{
    {
        auto r = bp::search_all("", bp::lit("XYZ"), bp::ws);
        int count = 0;
        for (auto subrange : r) {
            (void)subrange;
            ++count;
        }
        BOOST_TEST(count == 0);
    }
    {
        char const str[] = "aaXYZb";
        auto r = bp::search_all(str, bp::lit("XYZ"), bp::ws);
        int count = 0;
        int const offsets[] = {2, 5};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 1);
    }
    {
        char const str[] = "aaXYZbaabaXYZ";
        auto r = str | bp::search_all(bp::lit("XYZ"), bp::ws, bp::trace::off);
        int count = 0;
        int const offsets[] = {2, 5, 10, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str[] = "aaXYZbaabaXYZ";
        auto r = str | bp::search_all(bp::lit("XYZ"), bp::trace::off);
        int count = 0;
        int const offsets[] = {2, 5, 10, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str[] = "aaXYZbaabaXYZ";
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {2, 5, 10, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str[] = "aaXYZbaabaXYZXYZ";
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {2, 5, 10, 13, 13, 16};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str[] = "XYZaaXYZbaabaXYZXYZ";
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 3, 5, 8, 13, 16, 16, 19};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 3, 3, 6, 8, 11, 16, 19, 19, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
    {
        char const * str = "XYZXYZaaXYZbaabaXYZXYZ";
        auto r = bp::null_term(str) | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 3, 3, 6, 8, 11, 16, 19, 19, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
    {
        char const * str = "XYZXYZaaXYZbaabaXYZXYZ";
        auto const r = bp::null_term(str) | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 3, 3, 6, 8, 11, 16, 19, 19, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
}

// search_all_unicode
{
    {
        char const str_[] = "";
        auto str = str_ | bp::as_utf8;
        auto r = bp::search_all(str, bp::lit("XYZ"), bp::ws);
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
        auto r = bp::search_all(str, bp::lit("XYZ"), bp::ws);
        int count = 0;
        int const offsets[] = {2, 5};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 1);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf32;
        auto r = str | bp::search_all(bp::lit("XYZ"), bp::ws, bp::trace::off);
        int count = 0;
        int const offsets[] = {2, 5, 10, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf8;
        auto r = str | bp::search_all(bp::lit("XYZ"), bp::trace::off);
        int count = 0;
        int const offsets[] = {2, 5, 10, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {2, 5, 10, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str_[] = "aaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf32;
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {2, 5, 10, 13, 13, 16};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str_[] = "XYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf8;
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 3, 5, 8, 13, 16, 16, 19};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 3, 3, 6, 8, 11, 16, 19, 19, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
}

// doc_examples
{
    {
        namespace bp = boost::parser;
        auto result = bp::search("aaXYZq", bp::lit("XYZ"), bp::ws);
        assert(!result.empty());
        assert(
            std::string_view(result.begin(), result.end() - result.begin()) ==
            "XYZ");
        (void)result;
    }
    {
        auto r = "XYZaaXYZbaabaXYZXYZ" | bp::search_all(bp::lit("XYZ"));
        int count = 0;
        // Prints XYZ XYZ XYZ XYZ.
        for (auto subrange : r) {
            std::cout << std::string_view(subrange.begin(), subrange.end() - subrange.begin()) << " ";
            ++count;
        }
        std::cout << "\n";
        assert(count == 4);
    }
}

return boost::report_errors();
}
