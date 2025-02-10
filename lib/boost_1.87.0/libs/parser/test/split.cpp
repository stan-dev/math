/**
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/split.hpp>

#include <boost/core/lightweight_test.hpp>


namespace bp = boost::parser;

#if BOOST_PARSER_USE_CONCEPTS
namespace deduction {
    std::string str;
    auto const parser = bp::char_;
    auto const skip = bp::ws;

    auto deduced_1 = bp::split_view(str, parser, skip, bp::trace::on);
    auto deduced_2 = bp::split_view(str, parser, skip);
    auto deduced_3 = bp::split_view(str, parser, bp::trace::on);
    auto deduced_4 = bp::split_view(str, parser);
}
#endif

int main()
{

// split_
{
    {
        auto r = bp::split("", bp::lit("XYZ"), bp::ws);
        int count = 0;
        for (auto subrange : r) {
            (void)subrange;
            ++count;
        }
        BOOST_TEST(count == 0);
    }
    {
        char const str[] = "aaXYZb";
        auto r = bp::split(str, bp::lit("XYZ"), bp::ws);
        int count = 0;
        int const offsets[] = {0, 2, 5, 6};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str[] = "aaXYZbaabaXYZ";
        auto r = str | bp::split(bp::lit("XYZ"), bp::ws, bp::trace::off);
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str[] = "aaXYZbaabaXYZ";
        auto r = str | bp::split(bp::lit("XYZ"), bp::trace::off);
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str[] = "aaXYZbaabaXYZ";
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str[] = "aaXYZbaabaXYZXYZ";
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13, 16, 16};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str[] = "XYZaaXYZbaabaXYZXYZ";
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 0, 3, 5, 8, 13, 16, 16, 19, 19};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
    {
        char const str[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 0, 3, 3, 6, 8, 11, 16, 19, 19, 22, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 6);
    }
    {
        char const * str = "XYZXYZaaXYZbaabaXYZXYZ";
        auto r = bp::null_term(str) | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 0, 3, 3, 6, 8, 11, 16, 19, 19, 22, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 6);
    }
    {
        char const * str = "XYZXYZaaXYZbaabaXYZXYZ";
        auto const r = bp::null_term(str) | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 0, 3, 3, 6, 8, 11, 16, 19, 19, 22, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin() - str == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end() - str == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 6);
    }
}

// split_unicode
{
    {
        char const str_[] = "";
        auto str = str_ | bp::as_utf8;
        auto r = bp::split(str, bp::lit("XYZ"), bp::ws);
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
        auto r = bp::split(str, bp::lit("XYZ"), bp::ws);
        int count = 0;
        int const offsets[] = {0, 2, 5, 6};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 2);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf32;
        auto r = str | bp::split(bp::lit("XYZ"), bp::ws, bp::trace::off);
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf8;
        const auto r = str | bp::split(bp::lit("XYZ"), bp::trace::off);
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str_[] = "aaXYZbaabaXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 3);
    }
    {
        char const str_[] = "aaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf32;
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 2, 5, 10, 13, 13, 16, 16};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 4);
    }
    {
        char const str_[] = "XYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf8;
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 0, 3, 5, 8, 13, 16, 16, 19, 19};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 5);
    }
    {
        char const str_[] = "XYZXYZaaXYZbaabaXYZXYZ";
        auto str = str_ | bp::as_utf16;
        auto r = str | bp::split(bp::lit("XYZ"));
        int count = 0;
        int const offsets[] = {0, 0, 3, 3, 6, 8, 11, 16, 19, 19, 22, 22};
        for (auto subrange : r) {
            BOOST_TEST(subrange.begin().base() - str_ == offsets[count * 2 + 0]);
            BOOST_TEST(subrange.end().base() - str_ == offsets[count * 2 + 1]);
            ++count;
        }
        BOOST_TEST(count == 6);
    }
}

// doc_examples
{
    {
        auto r = "XYZaaXYZbaabaXYZXYZ" | bp::split(bp::lit("XYZ"));
        int count = 0;
        // Prints '' 'aa' 'baaba' '' ''.
        for (auto subrange : r) {
            std::cout << "'" << std::string_view(subrange.begin(), subrange.end() - subrange.begin()) << "' ";
            ++count;
        }
        std::cout << "\n";
        assert(count == 5);
    }
}

return boost::report_errors();
}
