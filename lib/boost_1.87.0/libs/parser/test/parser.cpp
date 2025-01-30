/**
 *   Copyright (C) 2018 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/parser.hpp>
#include <boost/parser/transcode_view.hpp>

#if __has_include(<boost/optional/optional.hpp>)
#define TEST_BOOST_OPTIONAL 1
#include <boost/optional/optional.hpp>
#else
#define TEST_BOOST_OPTIONAL 0
#endif

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;

#if TEST_BOOST_OPTIONAL
template<typename T>
constexpr bool boost::parser::enable_optional<boost::optional<T>> = true;
#endif


void msvc_only()
{

#if BOOST_PARSER_DETAIL_TEXT_USE_CONCEPTS && defined(_MSC_VER)
#include <algorithm>

    // msvc string_view
    {
        std::string_view sv = "text";
        std::wstring_view wsv = L"text";

        {
            auto r = detail::text::as_utf8(sv);
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = detail::text::as_utf8(wsv);
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
        {
            auto r = sv | detail::text::as_utf8;
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = wsv | detail::text::as_utf8;
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }

        {
            auto r = detail::text::as_utf16(sv);
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = detail::text::as_utf16(wsv);
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
        {
            auto r = sv | detail::text::as_utf16;
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = wsv | detail::text::as_utf16;
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }

        {
            auto r = detail::text::as_utf32(sv);
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = detail::text::as_utf32(wsv);
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
        {
            auto r = sv | detail::text::as_utf32;
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = wsv | detail::text::as_utf32;
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
    }

    // msvc string
    {
        std::string sv = "text";
        std::wstring wsv = L"text";

        {
            auto r = detail::text::as_utf8(sv);
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = detail::text::as_utf8(wsv);
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
        {
            auto r = sv | detail::text::as_utf8;
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = wsv | detail::text::as_utf8;
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }

        {
            auto r = detail::text::as_utf16(sv);
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = detail::text::as_utf16(wsv);
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
        {
            auto r = sv | detail::text::as_utf16;
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = wsv | detail::text::as_utf16;
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }

        {
            auto r = detail::text::as_utf32(sv);
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = detail::text::as_utf32(wsv);
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
        {
            auto r = sv | detail::text::as_utf32;
            std::string str = "    ";
            std::ranges::copy(r, str.begin());
            BOOST_TEST(str == sv);
        }
        {
            auto r = wsv | detail::text::as_utf32;
            std::wstring wstr = L"    ";
            std::ranges::copy(r, wstr.begin());
            BOOST_TEST(wstr == wsv);
        }
    }

#endif
}

namespace rule_construction_example {
    struct type_t
    {
        type_t() = default;
        explicit type_t(double x) : x_(x) {}
        // etc.

        double x_;
    };

    namespace bp = boost::parser;

    auto doubles_to_type = [](auto & ctx) {
        _val(ctx) = type_t(
            bp::get(_attr(ctx), bp::llong<0>{}) *
            bp::get(_attr(ctx), bp::llong<1>{}));
    };

    bp::rule<struct type_tag, type_t> type = "type";
    auto const type_def = (bp::double_ >> bp::double_)[doubles_to_type];
    BOOST_PARSER_DEFINE_RULES(type);
}

void rule_example()
{
    namespace bp = boost::parser;
    using namespace rule_construction_example;
    BOOST_TEST(bp::parse("3 4", type, bp::ws));
}

int main()
{

    msvc_only();
    rule_example();

    // basic
    {
        constexpr auto parser_1 = char_ >> char_;
        constexpr auto parser_2 = char_ >> char_ >> char_;
        constexpr auto parser_3 = char_ | char_;
        constexpr auto parser_4 = char_('a') | char_('b') | char_('c');
        constexpr auto parser_5 = char_('a') | char_('b') | eps;

        {
            std::string str = "a";
            BOOST_TEST(parse(str, char_));
            BOOST_TEST(!parse(str, char_('b')));
        }
        {
            std::string str = "a";
            char c = '\0';
            BOOST_TEST(parse(str, char_, c));
            BOOST_TEST(c == 'a');
            BOOST_TEST(!parse(str, char_('b')));
        }
        {
            std::string str = "b";
            char c = '\0';
            BOOST_TEST(parse(str, char_("ab"), c));
            BOOST_TEST(c == 'b');
            BOOST_TEST(!parse(str, char_("cd")));
        }
        {
            std::string str = "b";
            char c = '\0';
            std::string const pattern_1 = "ab";
            std::string const pattern_2 = "cd";
            BOOST_TEST(parse(str, char_(pattern_1), c));
            BOOST_TEST(c == 'b');
            BOOST_TEST(!parse(str, char_(pattern_2)));
        }
        {
            std::string str = "b";
            char c = '\0';
            BOOST_TEST(parse(str, char_('a', 'b'), c));
            BOOST_TEST(c == 'b');
            BOOST_TEST(!parse(str, char_('c', 'd')));
        }
        {
            std::string str = " ";
            BOOST_TEST(parse(str, blank));
            BOOST_TEST(!parse(str, lower));
        }
        {
            std::string str = "ab";
            BOOST_TEST(!parse(str, char_));
            {
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first, boost::parser::detail::text::null_sentinel, char_));
            }
            BOOST_TEST(parse(str, parser_1));
            BOOST_TEST(!parse(str, parser_2));
        }
        {
            std::string str = "ab";
            tuple<char, char> result;
            BOOST_TEST(parse(str, parser_1, result));
            using namespace boost::parser::literals;
            BOOST_TEST(get(result, 0_c) == 'b');
            BOOST_TEST(get(result, 1_c) == '\0');
        }
        {
            std::string str = "abc";
            BOOST_TEST(!parse(str, parser_1));
            {
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser_1));
            }
            BOOST_TEST(parse(str, parser_2));
        }
        {
            std::string str = "abc";
            tuple<char, char, char> result;
            BOOST_TEST(parse(str, parser_2, result));
            using namespace boost::parser::literals;
            BOOST_TEST(get(result, 0_c) == 'c');
            BOOST_TEST(get(result, 1_c) == '\0');
            BOOST_TEST(get(result, 2_c) == '\0');
        }
        {
            std::string str = "a";
            BOOST_TEST(parse(str, parser_3));
            BOOST_TEST(parse(str, parser_4));
        }
        {
            std::string str = "a";
            char c = '\0';
            BOOST_TEST(parse(str, parser_3, c));
            BOOST_TEST(c == 'a');
        }
        {
            std::string str = "a";
            char c = '\0';
            BOOST_TEST(parse(str, parser_4, c));
            BOOST_TEST(c == 'a');
        }
        {
            std::string str = "z";
            BOOST_TEST(parse(str, parser_3));
            BOOST_TEST(!parse(str, parser_4));
        }
        {
            std::string str = "a";
            BOOST_TEST(parse(str, parser_5));
        }
        {
            std::string str = "z";
            BOOST_TEST(!parse(str, parser_5));
            {
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser_5));
            }
        }
        {
            std::string str = "a";
            std::optional<char> c;
            BOOST_TEST(parse(str, parser_5, c));
            BOOST_TEST(c == 'a');
        }
        {
            std::string str = "z";
            std::optional<char> c;
            BOOST_TEST(!parse(str, parser_5, c));
        }
        {
            std::string str = "z";
            std::optional<char> c;
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                parser_5,
                c));
            BOOST_TEST(c == std::nullopt);
        }
#if TEST_BOOST_OPTIONAL
        {
            std::string str = "a";
            boost::optional<char> c;
            BOOST_TEST(parse(str, parser_5, c));
            BOOST_TEST(c == 'a');
        }
        {
            std::string str = "z";
            boost::optional<char> c;
            BOOST_TEST(!parse(str, parser_5, c));
        }
        {
            std::string str = "z";
            boost::optional<char> c;
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                parser_5,
                c));
            BOOST_TEST(c == boost::none);
        }
#endif
        {
            constexpr auto parser = *(char_ - "str") >> *string("str");
            std::string str = "somethingstrstrstr";
            auto result = boost::parser::parse(str, parser);
            BOOST_TEST(result);
            BOOST_TEST(
                *result ==
                std::vector<std::string>({"something", "str", "str", "str"}));
        }
    }

    // int_uint
    // clang-format off
    {
        {
            std::string str = "-42";
            int i = 0;
            BOOST_TEST(parse(str, int_, i));
            BOOST_TEST(i == -42);
        }
        // clang-format on
        {
            std::string str = "42";
            int i = 0;
            BOOST_TEST(parse(str, int_, i));
            BOOST_TEST(i == 42);
        }
        {
            std::string str = "-42";
            int i = 3;
            BOOST_TEST(!parse(str, uint_, i));
            BOOST_TEST(i == 0);
        }
        {
            std::string str = "42";
            int i = 0;
            BOOST_TEST(parse(str, uint_, i));
            BOOST_TEST(i == 42);
        }
        {
            std::string str = "-0042";
            int i = 0;
            parser_interface<int_parser<int, 10, 4, 4>> int4;
            BOOST_TEST(parse(str, int4, i));
            BOOST_TEST(i == -42);
        }
        {
            std::string str = "0042";
            unsigned int i = 0;
            parser_interface<uint_parser<unsigned int, 10, 4, 4>> uint4;
            BOOST_TEST(parse(str, uint4, i));
            BOOST_TEST(i == 42u);
        }
        }

        // bool_
        {{std::string str = "";
        bool b = false;
        BOOST_TEST(!parse(str, bool_, b));
        }
        {
            std::string str = "true";
            bool b = false;
            BOOST_TEST(parse(str, bool_, b));
            BOOST_TEST(b == true);
        }
        {
            std::string str = "false ";
            bool b = true;
            BOOST_TEST(!parse(str, bool_, b));
            BOOST_TEST(b == false);
        }
        {
            std::string str = "false ";
            bool b = true;
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first, boost::parser::detail::text::null_sentinel, bool_, b));
            BOOST_TEST(b == false);
        }
        {
            std::string str = "true ";
            auto r = boost::parser::detail::text::as_utf32(str);
            bool b = false;
            auto first = r.begin();
            auto const last = r.end();
            BOOST_TEST(prefix_parse(first, last, bool_, b));
            BOOST_TEST(b == true);
        }
        {
            std::string str = "false";
            auto r = boost::parser::detail::text::as_utf32(str);
            bool b = true;
            auto first = r.begin();
            auto const last = r.end();
            BOOST_TEST(prefix_parse(first, last, bool_, b));
            BOOST_TEST(b == false);
        }
        }

        // star
        {{constexpr auto parser = *char_;
        {
            std::string str = "";
            std::vector<char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<char>());
        }
        {
            std::string str = "a";
            std::vector<char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<char>({'a'}));
        }
        {
            std::string str = "ba";
            std::vector<char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'a'}));
        }
        }

        {
            constexpr auto parser = *char_('b');
            {
                std::string str = "";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>());
            }
            {
                std::string str = "b";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'b'}));
            }
            {
                std::string str = "bb";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'b', 'b'}));
            }
        }
        }

        // plus
        {{constexpr auto parser = +char_;

        {
            std::string str = "";
            std::vector<char> chars;
            BOOST_TEST(!parse(str, parser, chars));
        }
        {
            std::string str = "a";
            std::vector<char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<char>({'a'}));
        }
        {
            std::string str = "ba";
            std::vector<char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'a'}));
        }
        }

        {
            constexpr auto parser = +char_('b');

            {
                std::string str = "";
                std::vector<char> chars;
                BOOST_TEST(!parse(str, parser, chars));
            }
            {
                std::string str = "b";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'b'}));
            }
            {
                std::string str = "bb";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'b', 'b'}));
            }
        }
        }

        // star_and_plus_collapsing
        {{constexpr auto parser = +(+char_('b'));

        {
            std::string str = "";
            std::vector<char> chars;
            BOOST_TEST(!parse(str, parser, chars));
        }
        {
            std::string str = "b";
            std::vector<char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b'}));
        }
        {
            std::string str = "bb";
            std::vector<char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<char>({'b', 'b'}));
        }
        }

        {
            constexpr auto parser = **char_('z');

            {
                std::string str = "";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>());
            }
            {
                std::string str = "z";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z'}));
            }
            {
                std::string str = "zz";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
            }
        }

        {
            constexpr auto parser = +*char_('z');

            {
                std::string str = "";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>());
            }
            {
                std::string str = "z";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z'}));
            }
            {
                std::string str = "zz";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
            }
        }

        {
            constexpr auto parser = *+char_('z');

            {
                std::string str = "";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>());
            }
            {
                std::string str = "z";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z'}));
            }
            {
                std::string str = "zz";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
            }
        }
        }

        // action_
        {{std::string str = "";
        {
            std::stringstream ss;
            auto action = [&ss](auto & ctx) { ss << _attr(ctx); };
            auto parser = *char_('b')[action];
            BOOST_TEST(parse(str, parser));
            BOOST_TEST(ss.str() == "");
        }
        {
            str = "b";
            std::stringstream ss;
            auto action = [&ss](auto & ctx) { ss << _attr(ctx); };
            auto parser = *char_('b')[action];
            BOOST_TEST(parse(str, parser));
            BOOST_TEST(ss.str() == "b");
        }
        {
            str = "bb";
            std::stringstream ss;
            auto action = [&ss](auto & ctx) { ss << _attr(ctx); };
            auto parser = *char_('b')[action];
            BOOST_TEST(parse(str, parser));
            BOOST_TEST(parse(str, parser));
            BOOST_TEST(ss.str() == "bbbb");
        }
        }

        {
            {
                std::string str = "";
                std::stringstream ss;
                auto action = [&ss](auto & ctx) { ss << _attr(ctx); };
                auto parser = +char_('b')[action];
                BOOST_TEST(!parse(str, parser));
                BOOST_TEST(ss.str() == "");
            }
            {
                std::string str = "b";
                std::stringstream ss;
                auto action = [&ss](auto & ctx) { ss << _attr(ctx); };
                auto parser = +char_('b')[action];
                BOOST_TEST(parse(str, parser));
                BOOST_TEST(ss.str() == "b");
            }
            {
                std::string str = "bb";
                std::stringstream ss;
                auto action = [&ss](auto & ctx) { ss << _attr(ctx); };
                auto parser = +char_('b')[action];
                BOOST_TEST(parse(str, parser));
                BOOST_TEST(parse(str, parser));
                BOOST_TEST(ss.str() == "bbbb");
            }
        }
        }

        // star_as_string_or_vector
        {{constexpr auto parser = *char_('z');

        {
            std::string str = "";
            std::string chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == "");
        }
        {
            std::string str = "z";
            std::string chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == "z");
        }
        {
            std::string str = "zz";
            std::string chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == "zz");
        }
        }

        {
            constexpr auto parser = *char_('z');

            {
                std::string str = "";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>());
            }
            {
                std::string str = "z";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z'}));
            }
            {
                std::string str = "zz";
                std::vector<char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<char>({'z', 'z'}));
            }
        }

        {
            constexpr auto parser = *string("zs");

            {
                std::string str = "";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<std::string>{});

                {
                    std::optional<std::vector<std::string>> const chars =
                        parse(str, parser);
                    BOOST_TEST(chars);
                    BOOST_TEST(chars->empty());
                }
            }
            {
                std::string str = "z";
                {
                    std::vector<std::string> chars;
                    BOOST_TEST(!parse(str, parser, chars));
                    BOOST_TEST(chars == std::vector<std::string>{});
                }
                {
                    std::vector<std::string> chars;
                    auto first = str.c_str();
                    BOOST_TEST(prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser,
                        chars));
                    BOOST_TEST(chars == std::vector<std::string>{});
                }

                {
                    std::optional<std::vector<std::string>> const chars =
                        parse(str, parser);
                    BOOST_TEST(!chars);
                }
                {
                    auto first = str.c_str();
                    std::optional<std::vector<std::string>> const chars =
                        prefix_parse(
                            first,
                            boost::parser::detail::text::null_sentinel,
                            parser);
                    BOOST_TEST(chars);
                    BOOST_TEST(chars->empty());
                }
            }
            {
                std::string str = "zs";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<std::string>({"zs"}));

                {
                    std::optional<std::vector<std::string>> const chars =
                        parse(str, parser);
                    BOOST_TEST(chars);
                    BOOST_TEST(*chars == std::vector<std::string>({"zs"}));
                }
            }
            {
                std::string str = "zszs";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<std::string>({"zs", "zs"}));

                {
                    std::optional<std::vector<std::string>> const chars =
                        parse(str, parser);
                    BOOST_TEST(chars);
                    BOOST_TEST(
                        *chars == std::vector<std::string>({"zs", "zs"}));
                }
            }
        }

        {
            constexpr auto parser = *string("zs");

            {
                std::string str = "";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<std::string>());
            }
            {
                std::string str = "z";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, chars));
            }
            {
                std::string str = "z";
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    chars));
            }
            {
                std::string str = "zs";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<std::string>({"zs"}));
            }
            {
                std::string str = "zszs";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<std::string>({"zs", "zs"}));
            }
        }
        }

        // transform
        {
            int calls = 0;
            auto by_value_str_sum = [&](std::string s) {
                ++calls;
                std::transform(s.begin(), s.end(), s.begin(), [](auto ch) {
                    return ch - '0';
                });
                return std::accumulate(s.begin(), s.end(), 0);
            };
            auto cref_str_sum = [&](std::string const & s) {
                ++calls;
                int retval = 0;
                for (auto ch : s) {
                    retval += ch - '0';
                }
                return retval;
            };
            auto rv_ref_str_sum = [&](std::string && s) {
                ++calls;
                std::transform(s.begin(), s.end(), s.begin(), [](auto ch) {
                    return ch - '0';
                });
                return std::accumulate(s.begin(), s.end(), 0);
            };
            {
                constexpr auto parser = +char_;
                std::string str = "012345";
                {
                    auto result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == "012345");
                }
                {
                    calls = 0;
                    auto result =
                        parse(str, transform(by_value_str_sum)[parser]);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == 15);
                    BOOST_TEST(calls == 1);
                }
                {
                    calls = 0;
                    auto result = parse(str, transform(cref_str_sum)[parser]);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == 15);
                    BOOST_TEST(calls == 1);
                }
                {
                    calls = 0;
                    auto result = parse(str, transform(rv_ref_str_sum)[parser]);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == 15);
                    BOOST_TEST(calls == 1);
                }
            }
            {
                constexpr auto parser = +char_;
                std::string str = "012345";
                {
                    calls = 0;
                    auto result =
                        parse(str, omit[transform(by_value_str_sum)[parser]]);
                    BOOST_TEST(result);
                    BOOST_TEST(calls == 0);
                }
                {
                    calls = 0;
                    auto result =
                        parse(str, omit[transform(cref_str_sum)[parser]]);
                    BOOST_TEST(result);
                    BOOST_TEST(calls == 0);
                }
                {
                    calls = 0;
                    auto result =
                        parse(str, omit[transform(rv_ref_str_sum)[parser]]);
                    BOOST_TEST(result);
                    BOOST_TEST(calls == 0);
                }
            }
        }

        // transform_doc_example
        {
            //[ transform_directive_example
            auto str_sum = [&](std::string const & s) {
                int retval = 0;
                for (auto ch : s) {
                    retval += ch - '0';
                }
                return retval;
            };

            namespace bp = boost::parser;
            constexpr auto parser = +bp::char_;
            std::string str = "012345";

            auto result = bp::parse(str, bp::transform(str_sum)[parser]);
            assert(result);
            assert(*result == 15);
            static_assert(std::is_same_v<decltype(result), std::optional<int>>);
            //]
            (void)result;
        }

        // omit
        {{constexpr auto parser = omit[*+char_('z')];

        {
            std::string str = "";
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string str = "z";
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string str = "zz";
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string str = "";
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string str = "z";
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string str = "zz";
            BOOST_TEST(parse(str, parser));
        }
        }

        {
            constexpr auto parser = omit[*string("zs")];

            {
                std::string str = "";
                BOOST_TEST(parse(str, parser));
            }
            {
                std::string str = "z";
                BOOST_TEST(!parse(str, parser));
            }
            {
                std::string str = "z";
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first, boost::parser::detail::text::null_sentinel, parser));
            }
            {
                std::string str = "zs";
                BOOST_TEST(parse(str, parser));
            }
            {
                std::string str = "zszs";
                BOOST_TEST(parse(str, parser));
            }
        }
        }

        // repeat
        {{constexpr auto parser = repeat(2, 3)[string("zs")];

        {
            std::string str = "";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
        }
        {
            std::string str = "z";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
        }
        {
            std::string str = "zs";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
        }
        {
            std::string str = "zszs";
            std::vector<std::string> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"zs", "zs"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"zs", "zs"}));
            }
        }
        }
        {
            constexpr auto parser = *char_ >> eps >> *string("str");
            auto result = parse("abcdefg", parser);
            BOOST_TEST(result);
            BOOST_TEST(*result == std::vector<std::string>({"abcdefg"}));
        }
        {
            constexpr auto parser = *(char_ - "str") >> eps >> *string("str");
            auto result = parse("abcdefgstrstr", parser);
            BOOST_TEST(result);
            BOOST_TEST(
                *result == std::vector<std::string>({"abcdefg", "str", "str"}));
        }
        }

        // raw
        {
            {
                constexpr auto parser = raw[*string("zs")];
                using range_t = BOOST_PARSER_DETAIL_TEXT_SUBRANGE<
                    std::string::const_iterator>;

                {
                    std::string const str = "";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r.begin() == str.begin());
                    BOOST_TEST(r.end() == str.begin());
                }
                {
                    std::string const str = "z";
                    range_t r;
                    BOOST_TEST(!parse(str, parser, r));
                    BOOST_TEST(r.begin() == r.end());
                }
                {
                    std::string const str = "z";
                    range_t r;
                    auto first = str.begin();
                    BOOST_TEST(prefix_parse(first, str.end(), parser, r));
                    BOOST_TEST(r.begin() == str.begin());
                    BOOST_TEST(r.end() == str.begin());
                }
                {
                    std::string const str = "zs";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r.begin() == str.begin());
                    BOOST_TEST(r.end() == str.end());
                }
                {
                    std::string const str = "zszs";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r.begin() == str.begin());
                    BOOST_TEST(r.end() == str.end());
                }
                {
                    std::string const str = "";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(result->begin() == str.begin());
                    BOOST_TEST(result->end() == str.begin());
                }
                {
                    std::string const str = "z";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(!result);
                }
                {
                    std::string const str = "z";
                    auto first = str.begin();
                    std::optional<range_t> result =
                        prefix_parse(first, str.end(), parser);
                    BOOST_TEST(result);
                    BOOST_TEST(result->begin() == str.begin());
                    BOOST_TEST(result->end() == str.begin());
                }
                {
                    std::string const str = "zs";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(result->begin() == str.begin());
                    BOOST_TEST(result->end() == str.end());
                }
                {
                    std::string const str = "zszs";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(result->begin() == str.begin());
                    BOOST_TEST(result->end() == str.end());
                }
            }
        }

#if BOOST_PARSER_USE_CONCEPTS
        // string_view
        {
            {
                constexpr auto parser = string_view[*string("zs")];
                using range_t = std::string_view;

                {
                    std::string const str = "";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r == "");
                }
                {
                    std::string const str = "z";
                    range_t r;
                    BOOST_TEST(!parse(str, parser, r));
                    BOOST_TEST(r == "");
                }
                {
                    std::string const str = "z";
                    range_t r;
                    auto first = str.begin();
                    BOOST_TEST(prefix_parse(first, str.end(), parser, r));
                    BOOST_TEST(r == "");
                }
                {
                    std::string const str = "zs";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r == "zs");
                }
                {
                    std::string const str = "zszs";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r == "zszs");
                }
                {
                    std::string const str = "";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == "");
                }
                {
                    std::string const str = "z";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(!result);
                }
                {
                    std::string const str = "z";
                    auto first = str.begin();
                    std::optional<range_t> result =
                        prefix_parse(first, str.end(), parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == "");
                }
                {
                    std::string const str = "zs";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == "zs");
                }
                {
                    std::string const str = "zszs";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == "zszs");
                }
            }
            {
                constexpr auto parser = string_view[*string("zs")];
                using range_t = std::u32string_view;

                {
                    std::u32string const str = U"";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r == U"");
                }
                {
                    std::u32string const str = U"z";
                    range_t r;
                    BOOST_TEST(!parse(str, parser, r));
                    BOOST_TEST(r == U"");
                }
                {
                    std::u32string const str = U"z";
                    range_t r;
                    auto first = str.begin();
                    BOOST_TEST(prefix_parse(first, str.end(), parser, r));
                    BOOST_TEST(r == U"");
                }
                {
                    std::u32string const str = U"zs";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r == U"zs");
                }
                {
                    std::u32string const str = U"zszs";
                    range_t r;
                    BOOST_TEST(parse(str, parser, r));
                    BOOST_TEST(r == U"zszs");
                }
                {
                    std::u32string const str = U"";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == U"");
                }
                {
                    std::u32string const str = U"z";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(!result);
                }
                {
                    std::u32string const str = U"z";
                    auto first = str.begin();
                    std::optional<range_t> result =
                        prefix_parse(first, str.end(), parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == U"");
                }
                {
                    std::u32string const str = U"zs";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == U"zs");
                }
                {
                    std::u32string const str = U"zszs";
                    std::optional<range_t> result = parse(str, parser);
                    BOOST_TEST(result);
                    BOOST_TEST(*result == U"zszs");
                }
            }
        }
#endif

        // delimited
        {{constexpr auto parser = string("yay") % ',';

        {
            std::string str = "";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
        }
        {
            std::string str = "z";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
        }
        {
            std::string str = ",";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
        }
        {
            std::string str = ",yay";
            std::vector<std::string> chars;
            BOOST_TEST(!parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>{});

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
        }
        {
            std::string str = "yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == std::vector<std::string>({"yay"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            std::string str = "yayyay";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
            {
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            std::string str = "yay,";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(!chars);
            }
            {
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            std::string str = "yay,yay,yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(
                chars == std::vector<std::string>({"yay", "yay", "yay"}));

            {
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser);
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
        }
        }

        {
            constexpr auto parser = string("yay") % ',';
            {
                std::string str = "";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }

            {
                std::string str = "";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = "z";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }
            {
                std::string str = "z";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = ",";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }
            {
                std::string str = ",";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = " ,yay";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }
            {
                std::string str = " ,yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = ", yay";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }
            {
                std::string str = ", yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = ",yay ";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }
            {
                std::string str = ",yay ";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }

            {
                std::string str = " , yay ";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>{});
            }
            {
                std::string str = " , yay ";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = "yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }
            {
                std::string str = "yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
            }
            {
                std::string str = "yayyay";
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    char_(' '),
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }
            {
                std::string str = "yayyay";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
            }
            {
                std::string str = "yayyay";
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    char_(' '),
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }
            {
                std::string str = "yayyay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = "yayyay";
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser,
                        char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
            {
                std::string str = "yay,";
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
            }
            {
                std::string str = "yay,";
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    char_(' '),
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }
            {
                std::string str = "yay,";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = "yay,";
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser,
                        char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
            {
                std::string str = "yay,yay,yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay,yay,yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = " yay,yay,yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = " yay,yay,yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay ,yay,yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }

            {
                std::string str = "yay ,yay,yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay, yay,yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay, yay,yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay,yay ,yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }

            {
                std::string str = "yay,yay ,yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay,yay, yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }

            {
                std::string str = "yay,yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay,yay,yay ";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }

            {
                std::string str = "yay,yay,yay ";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = " yay , yay , yay ";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }

            {
                std::string str = " yay , yay , yay ";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
            {
                std::string str = "yay, yay, yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }

            {
                std::string str = "yay, yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(
                    *chars == std::vector<std::string>({"yay", "yay", "yay"}));
            }
        }

        {
            constexpr auto yay = string("yay") % ',';
            constexpr auto aww = string("aww") % ',';
            constexpr auto parser = yay >> ',' >> aww;

            {
                std::string str = "yay, yay, yay, aww, aww";
                auto result = parse(str, parser, ws);
                BOOST_TEST(result);
                BOOST_TEST(
                    *result ==
                    (tuple<std::vector<std::string>, std::vector<std::string>>(
                        std::vector<std::string>({"yay", "yay", "yay"}),
                        std::vector<std::string>({"aww", "aww"}))));
            }
            {
                std::string str = "yay, yay, yay , aww, aww";
                auto result = parse(str, parser, ws);
                BOOST_TEST(result);
                BOOST_TEST(
                    *result ==
                    (tuple<std::vector<std::string>, std::vector<std::string>>(
                        std::vector<std::string>({"yay", "yay", "yay"}),
                        std::vector<std::string>({"aww", "aww"}))));
            }
        }

        {
            constexpr auto yay = string("yay") % ',';
            constexpr auto aww = string("aww") % ',';
            constexpr auto parser = raw[yay] >> ',' >> raw[aww];

            using namespace boost::parser::literals;

            {
                std::string str = "yay, yay, yay, aww, aww";
                auto result = parse(str, parser, ws);
                BOOST_TEST(result);
                auto subrange_0 = boost::parser::get(*result, 0_c);
                auto subrange_1 = boost::parser::get(*result, 1_c);
                BOOST_TEST(subrange_0.begin() == str.begin());
                BOOST_TEST(subrange_0.end() == str.begin() + 13);
                BOOST_TEST(subrange_1.begin() == str.begin() + 15);
                BOOST_TEST(subrange_1.end() == str.begin() + 23);
            }
            {
                std::string str = "yay, yay, yay , aww, aww";
                auto result = parse(str, parser, ws);
                BOOST_TEST(result);
                auto subrange_0 = boost::parser::get(*result, 0_c);
                auto subrange_1 = boost::parser::get(*result, 1_c);
                BOOST_TEST(subrange_0.begin() == str.begin());
                BOOST_TEST(subrange_0.end() == str.begin() + 13);
                BOOST_TEST(subrange_1.begin() == str.begin() + 16);
                BOOST_TEST(subrange_1.end() == str.begin() + 24);
            }
        }

        {
            constexpr auto yay = *string("yay");
            constexpr auto aww = *string("aww");
            constexpr auto parser = raw[yay] >> ',' >> raw[aww];

            {
                std::string str = "yay yay yay, aww aww";
                auto result = parse(str, parser, ws);
                BOOST_TEST(result);
                auto subrange_0 = boost::parser::get(*result, llong<0>{});
                auto subrange_1 = boost::parser::get(*result, llong<1>{});
                BOOST_TEST(subrange_0.begin() == str.begin());
                BOOST_TEST(subrange_0.end() == str.begin() + 11);
                BOOST_TEST(subrange_1.begin() == str.begin() + 13);
                BOOST_TEST(subrange_1.end() == str.begin() + 20);
            }
            {
                std::string str = "yay yay yay , aww aww";
                auto result = parse(str, parser, ws);
                BOOST_TEST(result);
                auto subrange_0 = boost::parser::get(*result, llong<0>{});
                auto subrange_1 = boost::parser::get(*result, llong<1>{});
                BOOST_TEST(subrange_0.begin() == str.begin());
                BOOST_TEST(subrange_0.end() == str.begin() + 11);
                BOOST_TEST(subrange_1.begin() == str.begin() + 14);
                BOOST_TEST(subrange_1.end() == str.begin() + 21);
            }
        }
        }

        // lexeme
        {{constexpr auto parser = lexeme[string("yay") % ','];

        {
            std::string str = "yay, yay, yay";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    char_(' '),
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                std::string str = "yay, yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = "yay, yay, yay";
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser,
                        char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        {
            std::string str = " yay, yay, yay";
            {
                std::vector<std::string> chars;
                BOOST_TEST(!parse(str, parser, char_(' '), chars));
            }
            {
                std::vector<std::string> chars;
                auto first = str.c_str();
                BOOST_TEST(prefix_parse(
                    first,
                    boost::parser::detail::text::null_sentinel,
                    parser,
                    char_(' '),
                    chars));
                BOOST_TEST(chars == std::vector<std::string>({"yay"}));
            }

            {
                std::string str = " yay, yay, yay";
                std::optional<std::vector<std::string>> const chars =
                    parse(str, parser, char_(' '));
                BOOST_TEST(!chars);
            }
            {
                std::string str = " yay, yay, yay";
                auto first = str.c_str();
                std::optional<std::vector<std::string>> const chars =
                    prefix_parse(
                        first,
                        boost::parser::detail::text::null_sentinel,
                        parser,
                        char_(' '));
                BOOST_TEST(chars);
                BOOST_TEST(*chars == std::vector<std::string>({"yay"}));
            }
        }
        }

        {
            constexpr auto parser = lexeme[skip[string("yay") % ',']];

            {
                std::string str = "yay, yay, yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));

                {
                    std::string str = "yay, yay, yay";
                    std::optional<std::vector<std::string>> const chars =
                        parse(str, parser, char_(' '));
                    BOOST_TEST(chars);
                    BOOST_TEST(
                        *chars ==
                        std::vector<std::string>({"yay", "yay", "yay"}));
                }
            }
            {
                std::string str = " yay, yay, yay";
                std::vector<std::string> chars;
                BOOST_TEST(parse(str, parser, char_(' '), chars));
                BOOST_TEST(
                    chars == std::vector<std::string>({"yay", "yay", "yay"}));

                {
                    std::string str = " yay, yay, yay";
                    std::optional<std::vector<std::string>> const chars =
                        parse(str, parser, char_(' '));
                    BOOST_TEST(chars);
                    BOOST_TEST(
                        *chars ==
                        std::vector<std::string>({"yay", "yay", "yay"}));
                }
            }
        }

        {
            auto const parser = string("abc");

            // Follows the parser used in transform_replace().
            auto before = [&](auto & ctx) {};
            auto after = [&](auto & ctx) {};
            auto const search_parser =
                omit[*(char_ - parser)] >>
                -lexeme[eps[before] >> skip[parser] >> eps[after]];

            {
                std::string str = "abc";
                std::optional<std::string> result;
                BOOST_TEST(parse(str, search_parser, char_(' '), result));
                BOOST_TEST(*result == "abc");

                {
                    std::string str = "abc";
                    auto first = detail::text::detail::begin(str);
                    auto last = detail::text::detail::end(str);
                    auto const result =
                        prefix_parse(first, last, search_parser, char_(' '));
                    static_assert(
                        std::is_same_v<
                            decltype(result),
                            std::optional<std::optional<std::string>> const>);
                    BOOST_TEST(result);
                    BOOST_TEST(**result == "abc");
                }
            }
            {
                std::string str = " abc";
                std::optional<std::string> result;
                BOOST_TEST(parse(str, search_parser, char_(' '), result));
                BOOST_TEST(*result == "abc");

                {
                    std::string str = " abc";
                    auto const result = parse(str, search_parser, char_(' '));
                    static_assert(
                        std::is_same_v<
                            decltype(result),
                            std::optional<std::optional<std::string>> const>);
                    BOOST_TEST(result);
                    BOOST_TEST(**result == "abc");
                }
            }
        }

        {
            auto const parser = int_ % ',';

            // Follows the parser used in transform_replace().
            auto before = [&](auto & ctx) {};
            auto after = [&](auto & ctx) {};
            auto const search_parser =
                omit[*(char_ - parser)] >>
                -lexeme[eps[before] >> skip[parser] >> eps[after]];

            {
                std::string str = "1, 2, 4";
                std::optional<std::vector<int>> result;
                BOOST_TEST(parse(str, search_parser, char_(' '), result));
                BOOST_TEST(*result == std::vector<int>({1, 2, 4}));

                {
                    std::string str = "1, 2, 4";
                    auto const result = parse(str, search_parser, char_(' '));
                    static_assert(std::is_same_v<
                                  decltype(result),
                                  std::optional<
                                      std::optional<std::vector<int>>> const>);
                    BOOST_TEST(result);
                    BOOST_TEST(**result == std::vector<int>({1, 2, 4}));
                }
            }
            {
                std::string str = " 1, 2, 4";
                std::optional<std::vector<int>> result;
                BOOST_TEST(parse(str, search_parser, char_(' '), result));
                BOOST_TEST(*result == std::vector<int>({1, 2, 4}));

                {
                    std::string str = " 1, 2, 4";
                    auto const result = parse(str, search_parser, char_(' '));
                    static_assert(std::is_same_v<
                                  decltype(result),
                                  std::optional<
                                      std::optional<std::vector<int>>> const>);
                    BOOST_TEST(result);
                    BOOST_TEST(**result == std::vector<int>({1, 2, 4}));
                }
            }
        }
        }

        // skip
        {{constexpr auto parser = skip(char_(' '))[string("yay") % ','];

        {
            std::string str = "yay, yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(
                chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            std::string str = "yay, yay, yay";
            std::optional<std::vector<std::string>> const chars =
                parse(str, parser);
            BOOST_TEST(chars);
            BOOST_TEST(
                *chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            std::string str = " yay, yay, yay";
            std::vector<std::string> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(
                chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        {
            std::string str = " yay, yay, yay";
            std::optional<std::vector<std::string>> const chars =
                parse(str, parser);
            BOOST_TEST(chars);
            BOOST_TEST(
                *chars == std::vector<std::string>({"yay", "yay", "yay"}));
        }
        }
        }

        // combined_seq_and_or
        {{constexpr auto parser = char_('a') >> char_('b') >> char_('c') |
                                  char_('x') >> char_('y') >> char_('z');
        using tup = tuple<char, char, char>;

        {
            std::string str = "abc";
            tuple<char, char, char> chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == tup('c', '\0', '\0'));
        }

        {
            std::string str = "abc";
            std::optional<std::string> const chars = parse(str, parser);
            BOOST_TEST(chars);
            BOOST_TEST(*chars == "abc");
        }

        {
            std::string str = "xyz";
            tup chars;
            BOOST_TEST(parse(str, parser, chars));
            BOOST_TEST(chars == tup('z', '\0', '\0'));
        }
        }

        {
            constexpr auto parser = char_('a') >> string("b") >> char_('c') |
                                    char_('x') >> string("y") >> char_('z');
            {
                std::string str = "abc";
                std::string chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == "abc");
            }

            {
                std::string str = "abc";
                std::optional<std::string> const chars = parse(str, parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == "abc");
            }

            {
                std::string str = "xyz";
                std::string chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == "xyz");
            }
        }

        {
            constexpr auto parser = char_('a') >> char_('b') >> char_('c') |
                                    char_('x') >> char_('y') >> char_('z');
            using tup = tuple<char, char, char>;

            {
                std::string str = "abc";
                tuple<std::any, std::any, std::any> chars;
                BOOST_TEST(parse(str, parser, chars));
            }

            {
                std::string str = "xyz";
                tuple<char, char, char> chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == tup('z', '\0', '\0'));
            }
        }

        {
            constexpr auto parser = !char_('a');
            std::string str = "a";
            BOOST_TEST(!parse(str, parser));
        }

        {
            constexpr auto parser = &char_('a');
            std::string str = "a";
            BOOST_TEST(!parse(str, parser));
        }
        {
            constexpr auto parser = &char_('a');
            std::string str = "a";
            auto first = str.c_str();
            BOOST_TEST(prefix_parse(
                first, boost::parser::detail::text::null_sentinel, parser));
        }

        {
#if defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-shift-op-parentheses"
#endif
            constexpr auto parser = (char_('a') >> string("b") > char_('c')) |
                                    (char_('x') >> string("y") >> char_('z'));
#if defined(__clang__)
#pragma GCC diagnostic pop
#endif
            {
                std::string str = "abc";
                std::string chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == "abc");
            }

            {
                std::string str = "abc";
                std::optional<std::string> const chars = parse(str, parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == "abc");
            }

            {
                std::string str = "xyz";
                std::string chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == "xyz");
            }

            {
                std::string str = "abz";
                std::string chars;
                rethrow_error_handler eh;
                BOOST_TEST_THROWS(
                    parse(str, with_error_handler(parser, eh), chars),
                    parse_error<std::string::const_iterator>);
            }

            {
                std::string str = "abz";
                std::string chars;
                BOOST_TEST(!parse(str, parser, chars));
            }

            {
                std::string str = "abz";
                std::string chars;
                stream_error_handler eh("simple_parser.cpp");
                BOOST_TEST(!parse(str, with_error_handler(parser, eh), chars));
            }

            {
                std::string str = "ab";
                std::string chars;
                stream_error_handler eh("simple_parser.cpp");
                BOOST_TEST(!parse(str, with_error_handler(parser, eh), chars));
            }
        }

        {
#if defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-shift-op-parentheses"
#endif
            constexpr auto parser = (char_('a') >> string("b") > char_('c')) |
                                    (char_('x') >> string("y") >> char_('z'));
#if defined(__clang__)
#pragma GCC diagnostic pop
#endif
            {
                std::string str = "abc";
                std::string chars;
                BOOST_TEST(parse(str, parser, chars));
                BOOST_TEST(chars == "abc");
            }

            {
                std::string str = "abc";
                std::optional<std::string> const chars = parse(str, parser);
                BOOST_TEST(chars);
                BOOST_TEST(*chars == "abc");
            }

            {
                std::string str = "xyz";
                std::string chars;
                BOOST_TEST(parse(str, parser, chars, trace::on));
                BOOST_TEST(chars == "xyz");
            }
        }

        {
            constexpr auto parser = int_ >> -(lit('a') | 'b');
            auto result = parse("34b", parser);
            BOOST_TEST(result);
            BOOST_TEST(*result == 34);
        }
        }

        // eol_
        {{constexpr auto parser = eol;

        {
            std::string str = "y";
            BOOST_TEST(!parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000d\u000a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000b";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000c";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000d";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u0085";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2028";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2029";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        }
        }

        // ws_
        {{constexpr auto parser = ws;

        {
            std::string str = "y";
            BOOST_TEST(!parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u0009";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000d\u000a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000b";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000c";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u000d";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u0085";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u00a0";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u1680";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2000";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2001";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2002";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2003";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2004";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2005";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2006";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2007";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2008";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2009";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u200a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2028";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u2029";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u202F";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u205F";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u3000";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser));
        }
        }
        }

        // blank_
        {{constexpr auto parser = blank;
        constexpr auto alt_parser = ws - eol;

        {
            std::string str = "y";
            BOOST_TEST(!parse(str, parser));
        }
        {
            std::string s = (char const *)u8"\u0009";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u000a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u000d\u000a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u000b";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u000c";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u000d";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u0085";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u00a0";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u1680";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2000";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2001";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2002";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2003";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2004";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2005";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2006";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2007";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2008";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2009";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u200a";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2028";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u2029";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u202F";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u205F";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        {
            std::string s = (char const *)u8"\u3000";
            auto str = boost::parser::detail::text::as_utf8(s);
            BOOST_TEST(parse(str, parser) == parse(str, alt_parser));
        }
        }
        }

        // digit_
        {
            constexpr auto parser = +digit;

            std::u32string str = U"a/09:\x0659\x0660\x0669\x066a";
            std::vector<uint32_t> result;
            BOOST_TEST(parse(str, parser, char_ - digit, result));
            BOOST_TEST(
                result == std::vector<uint32_t>({'0', '9', 0x0660, 0x0669}));
        }

        // hex_digit_
        {
            constexpr auto parser = +hex_digit;

            std::u32string str = U"a/09:A\uff0f\uff10\uff19\uff1a";
            std::vector<uint32_t> result;
            BOOST_TEST(parse(str, parser, char_ - hex_digit, result));
            BOOST_TEST(
                result == std::vector<uint32_t>(
                              {'a', '0', '9', 'A', U'\uff10', U'\uff19'}));
        }

        // control_
        {
            constexpr auto parser = +control;

            std::u32string str = U"\u0001\u001f\u0020\u007e\u007f\u009f\u00a0";
            std::vector<uint32_t> result;
            BOOST_TEST(parse(str, parser, char_ - control, result));
            BOOST_TEST(
                result == std::vector<uint32_t>({1, 0x001f, 0x007f, 0x009f}));
        }

        // punct_
        {
            auto parser = +punct;

            std::u32string str = U"\u0020\u0021\u0fda\u0fdb";
            std::vector<uint32_t> result;
            BOOST_TEST(parse(str, parser, char_ - punct, result));
            BOOST_TEST(result == std::vector<uint32_t>({0x21, 0xfda}));
        }

        // lower_
        {
            auto parser = +lower;

            std::u32string str = U"aA\u016F\u0170";
            std::vector<uint32_t> result;
            BOOST_TEST(parse(str, parser, char_ - lower, result));
            BOOST_TEST(result == std::vector<uint32_t>({'a', 0x16f}));
        }

        // upper_
        {
            auto parser = +upper;

            std::u32string str = U"aA\u0105\u0106";
            std::vector<uint32_t> result;
            BOOST_TEST(parse(str, parser, char_ - upper, result));
            BOOST_TEST(result == std::vector<uint32_t>({'A', 0x106}));
        }

        // no_need_for_sprit_2_hold_directive
        {
            namespace bp = boost::parser;

            std::vector<int> v;
            auto result = bp::parse(
                "1 2",
                bp::repeat(3)[bp::int_] | repeat(2)[bp::int_] >> bp::attr(0),
                bp::ws,
                v,
                bp::trace::on);
            BOOST_TEST(result);

            BOOST_TEST(v.size() == 3u);
            BOOST_TEST(v == std::vector<int>({1, 2, 0}));
        }

        // raw_doc_example
        {
            namespace bp = boost::parser;
            auto int_parser =
                bp::int_ % ','; // ATTR(int_parser) is std::vector<int>
            auto subrange_parser =
                bp::raw[int_parser]; // ATTR(subrange_parser) is a subrange

            // Parse using int_parser, generating integers.
            auto ints = bp::parse("1, 2, 3, 4", int_parser, bp::ws);
            assert(ints);
            assert(*ints == std::vector<int>({1, 2, 3, 4}));

            // Parse again using int_parser, but this time generating only the
            // subrange matched by int_parser.  (prefix_parse() allows matches
            // that don't consume the entire input.)
            auto const str = std::string("1, 2, 3, 4, a, b, c");
            auto first = str.begin();
            auto range =
                bp::prefix_parse(first, str.end(), subrange_parser, bp::ws);
            assert(range);
            assert(range->begin() == str.begin());
            assert(range->end() == str.begin() + 10);

            static_assert(
                std::is_same_v<
                    decltype(range),
                    std::optional<
                        BOOST_PARSER_SUBRANGE<std::string::const_iterator>>>);

#if defined(__cpp_char8_t)
            auto const u8str = std::u8string(u8"1, 2, 3, 4, a, b, c");
            auto u8first = u8str.begin();
            auto u8range =
                bp::prefix_parse(u8first, u8str.end(), subrange_parser, bp::ws);
            assert(u8range);
            assert(u8range->begin().base() == u8str.begin());
            assert(u8range->end().base() == u8str.begin() + 10);
#endif
        }


#if BOOST_PARSER_USE_CONCEPTS
        // string_view_doc_example
        {
            namespace bp = boost::parser;
            auto int_parser =
                bp::int_ % ','; // ATTR(int_parser) is std::vector<int>
            auto sv_parser = bp::string_view[int_parser]; // ATTR(subrange_parser)
                                                          // is a string_view

            auto const str = std::string("1, 2, 3, 4, a, b, c");
            auto first = str.begin();
            auto sv1 = bp::prefix_parse(first, str.end(), sv_parser, bp::ws);
            assert(sv1);
            assert(*sv1 == str.substr(0, 10));

            static_assert(
                std::is_same_v<decltype(sv1), std::optional<std::string_view>>);

            auto sv2 =
                bp::parse("1, 2, 3, 4" | bp::as_utf32, sv_parser, bp::ws);
            assert(sv2);
            assert(*sv2 == "1, 2, 3, 4");

            static_assert(
                std::is_same_v<decltype(sv2), std::optional<std::string_view>>);
        }
#endif

        // variant_compat_example
        {
            struct key_value
            {
                int key;
                double value;
            };

            namespace bp = boost::parser;
            key_value kv;
            bp::parse("42 13.0", bp::int_ >> bp::double_, kv); // Ok.
#if 0
            std::variant<key_value, double> kv_or_d;
            bp::parse("42 13.0", bp::int_ >> bp::double_, kv_or_d); // Error: ill-formed!
#endif
        }


        // utf_iterator_copy_ctor
        {
            namespace bp = boost::parser;

            char mut_chars[] = "foo";
            char const const_chars[] = "bar";

            auto mut_view = mut_chars | bp::as_utf8;
            auto const_view = const_chars | bp::as_utf8;

            auto mut_it = mut_view.begin();
            auto mut_last = mut_view.end();
            auto const_it = const_view.begin();
            auto const_last = const_view.begin();

            const_it = mut_it;
            const_last = mut_last;

            std::string copy;
            copy.resize(3);
            std::copy(const_it, const_last, copy.begin());
            BOOST_TEST(copy == "foo");
        }

        // detail_printing
        {
            // print_printable(char > 127)
            {
                std::ostringstream oss;
                detail::print_printable(oss, char(130));
                BOOST_TEST(oss.str() == "'\\x82'");
            }

            // print_printable(char32_t)
            {
                std::ostringstream oss;
                detail::print_printable(oss, U'a');
                BOOST_TEST(oss.str() == "U'a'");
            }
            {
                std::ostringstream oss;
                detail::print_printable(oss, U'ß');
                BOOST_TEST(oss.str() == "U'\\xdf'");
            }

            // print(tuple)
            {
                std::ostringstream oss;
                tuple<int, float> tup(42, 13.8);
                detail::print(oss, tup);
                BOOST_TEST(oss.str() == "(42, 13.8)");
            }

            // print(optional)
            {
                std::ostringstream oss;
                detail::print(oss, std::optional<int>());
                BOOST_TEST(oss.str() == "<<empty>>");
            }
            {
                std::ostringstream oss;
                detail::print(oss, std::optional(42));
                BOOST_TEST(oss.str() == "42");
            }

            // print(variant)
            {
                std::ostringstream oss;
                detail::print(oss, std::variant<int, double>());
                BOOST_TEST(oss.str() == "<<variant>>");
            }
        }

        return boost::report_errors();
        }
