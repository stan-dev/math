// Copyright (C) 2018 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>


using namespace boost::parser;

void compile_attribute_non_unicode()
{
    // range
    {
        char const r[] = "";

        {
            constexpr auto parser = char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<char>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          char>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
    }
    // pointer-as-range (covers iter/sent case, as that's what gets used
    // internally)
    {
        char const * r = "";

        {
            constexpr auto parser = char_;
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<char>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          char>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
    }
    // iter/iter
    {
        char const chars[] = "";
        auto first = std::begin(chars);
        auto const last = std::end(chars);

        {
            constexpr auto parser = char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<char>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          char>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
    }
}

void compile_attribute_unicode_utf8()
{
#if defined(__cpp_char8_t)
    {
        char8_t const chars[] = u8"";
        auto first = std::begin(chars);
        auto const last = std::end(chars);

        {
            constexpr auto parser = char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          char32_t>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
    }
#endif
    {
        char const chars[] = "";
        auto const r = boost::parser::detail::text::as_utf8(chars);

        {
            constexpr auto parser = char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          char32_t>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
    }
#if defined(__cpp_char8_t)
    {
        char8_t const r[] = u8"";

        {
            constexpr auto parser = char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          char32_t>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(parse(r, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(r), decltype(parser)>,
                          std::string>);
        }
    }
    {
        char8_t const * r = u8"";

        {
            constexpr auto parser = char_;
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(null_term(r)), decltype(parser)>,
                          char32_t>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(null_term(r)), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(null_term(r)), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(null_term(r)), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(parse(null_term(r), parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(null_term(r)), decltype(parser)>,
                          std::string>);
        }
    }
#endif
}

void compile_attribute_unicode_utf32()
{
    {
        char32_t const chars[] = {0};
        auto first = std::begin(chars);
        auto const last = std::end(chars);

        {
            constexpr auto parser = char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          char32_t>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(prefix_parse(first, last, parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<
                              decltype(BOOST_PARSER_SUBRANGE(first, last)),
                              decltype(parser)>,
                          std::string>);
        }
    }
    {
        {
            constexpr auto parser = char_;
            using attr_t = decltype(parse(U"", parser));
            static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(U""), decltype(parser)>,
                          char32_t>);
        }
        {
            constexpr auto parser = *char_;
            using attr_t = decltype(parse(U"", parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(U""), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo");
            using attr_t = decltype(parse(U"", parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(U""), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = char_ >> string("foo");
            using attr_t = decltype(parse(U"", parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(U""), decltype(parser)>,
                          std::string>);
        }
        {
            constexpr auto parser = string("foo") >> char_;
            using attr_t = decltype(parse(U"", parser));
            static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
            static_assert(std::is_same_v<
                          attribute_t<decltype(U""), decltype(parser)>,
                          std::string>);
        }
    }
}

[[maybe_unused]] rule<class test_rule, std::string> const test_rule =
    "test_rule";
auto const test_rule_def = +char_;
BOOST_PARSER_DEFINE_RULES(test_rule);

[[maybe_unused]] rule<class ints, std::vector<int>> const ints = "ints";
auto twenty_zeros = [](auto & ctx) { _val(ctx).resize(20, 0); };
auto push_back = [](auto & ctx) { _val(ctx).push_back(_attr(ctx)); };
auto const ints_def = lit("20-zeros")[twenty_zeros] | +int_[push_back];
BOOST_PARSER_DEFINE_RULES(ints);

void compile_attribute_sentinel()
{
    char const * chars = "";
    auto first = chars;
    boost::parser::detail::text::null_sentinel_t last;

    {
        constexpr auto parser = eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = eol;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = eoi;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = attr(3.0);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<double>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      double>);
    }
    {
        constexpr auto parser = attr('c');
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        constexpr auto parser = cp;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char32_t>);
    }
    {
        constexpr auto parser = cu;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = 'c'_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8'c'_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U'c'_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = "str"_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8"str"_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U"str"_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = 'c'_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8'c'_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U'c'_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = "str"_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8"str"_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U"str"_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }
    {
        constexpr auto parser = lit('c');
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = lit(U'c');
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = lit("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        std::string str = "str";
        auto parser = lit(str);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = bool_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<bool>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      bool>);
    }
    {
        constexpr auto parser = bin;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = oct;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = hex;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = ushort_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned short>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned short>);
    }
    {
        constexpr auto parser = uint_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = ulong_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned long>);
    }
    {
        constexpr auto parser = ulong_long;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<unsigned long long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned long long>);
    }
    {
        constexpr auto parser = short_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<short>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      short>);
    }
    {
        constexpr auto parser = int_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      int>);
    }
    {
        constexpr auto parser = long_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      long>);
    }
    {
        constexpr auto parser = long_long;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<long long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      long long>);
    }
    {
        constexpr auto parser = float_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<float>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      float>);
    }
    {
        constexpr auto parser = double_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<double>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      double>);
    }
    {
        symbols<float> parser;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<float>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      float>);
    }
    {
        auto const c = [](auto & ctx) { return true; };
        constexpr auto parser = if_(c)[double_];
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<double>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      double>);
    }
    {
        auto const x = [](auto & ctx) { return 2; };
        auto parser = switch_(x)(0, double_)(2, int_);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::variant<double, int>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::variant<double, int>>);
    }

    {
        constexpr auto parser = int_ | int_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      int>);
    }
    {
        constexpr auto parser = double_ | int_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::variant<double, int>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::variant<double, int>>);
    }
    {
        constexpr auto parser = double_ | int_ | eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<std::optional<std::variant<double, int>>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::optional<std::variant<double, int>>>);
    }

    {
        constexpr auto parser = *cu >> string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::string, std::string>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      tuple<std::string, std::string>>);
    }
    {
        constexpr auto parser = cu >> string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }

    {
        auto a = [](auto & ctx) {};
        constexpr auto parser = cu[a];
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }

    {
        using attr_t = decltype(prefix_parse(first, last, test_rule));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
#if 0 // Intentionally ill-formed
        std::vector<char> vec;
        auto result = prefix_parse(first, last, test_rule, vec);
#endif
    }

    {
        using attr_t = decltype(prefix_parse(first, last, ints));
        static_assert(std::is_same_v<attr_t, std::optional<std::vector<int>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(ints)>,
                      std::vector<int>>);
    }
}

void compile_attribute()
{
    compile_attribute_non_unicode();
    compile_attribute_unicode_utf8();
    compile_attribute_unicode_utf32();
    compile_attribute_sentinel();

    char const chars[] = "";
    auto first = std::begin(chars);
    auto const last = std::end(chars);

    {
        constexpr auto parser = eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = eol;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = eoi;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = attr(3.0);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<double>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      double>);
    }
    {
        constexpr auto parser = attr('c');
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        constexpr auto parser = cp;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char32_t>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char32_t>);
    }
    {
        constexpr auto parser = cu;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = 'c'_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8'c'_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U'c'_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = "str"_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8"str"_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U"str"_l;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = 'c'_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8'c'_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U'c'_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<char>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      char>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = "str"_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = u8"str"_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }
    {
        using namespace boost::parser::literals;
        constexpr auto parser = U"str"_p;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }
    {
        constexpr auto parser = lit('c');
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = lit(U'c');
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = lit("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        std::string str = "str";
        auto parser = lit(str);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }
    {
        constexpr auto parser = bool_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<bool>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      bool>);
    }
    {
        constexpr auto parser = bin;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = oct;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = hex;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = ushort_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned short>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned short>);
    }
    {
        constexpr auto parser = uint_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned int>);
    }
    {
        constexpr auto parser = ulong_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<unsigned long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned long>);
    }
    {
        constexpr auto parser = ulong_long;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<unsigned long long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      unsigned long long>);
    }
    {
        constexpr auto parser = short_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<short>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      short>);
    }
    {
        constexpr auto parser = int_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      int>);
    }
    {
        constexpr auto parser = long_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      long>);
    }
    {
        constexpr auto parser = long_long;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<long long>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      long long>);
    }
    {
        constexpr auto parser = float_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<float>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      float>);
    }
    {
        constexpr auto parser = double_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<double>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      double>);
    }
    {
        symbols<float> parser;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<float>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      float>);
    }
    {
        auto const c = [](auto & ctx) { return true; };
        constexpr auto parser = if_(c)[double_];
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<double>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      double>);
    }
    {
        auto const x = [](auto & ctx) { return 2; };
        auto parser = switch_(x)(0, double_)(2, int_);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::variant<double, int>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::variant<double, int>>);
    }

    {
        constexpr auto parser = int_ | int_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<int>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      int>);
    }
    {
        constexpr auto parser = double_ | int_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::variant<double, int>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::variant<double, int>>);
    }
    {
        constexpr auto parser = double_ | int_ | eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<std::optional<std::variant<double, int>>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::optional<std::variant<double, int>>>);
    }

    {
        constexpr auto parser = *cu >> string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::string, std::string>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      tuple<std::string, std::string>>);
    }
    {
        constexpr auto parser = cu >> string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }

    {
        auto a = [](auto & ctx) {};
        constexpr auto parser = cu[a];
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, bool>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      none>);
    }

    {
        using attr_t = decltype(prefix_parse(first, last, test_rule));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
#if 0 // Intentionally ill-formed
        std::vector<char> vec;
        auto result = prefix_parse(first, last, test_rule, vec);
#endif
    }

    {
        using attr_t = decltype(prefix_parse(first, last, ints));
        static_assert(std::is_same_v<attr_t, std::optional<std::vector<int>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(ints)>,
                   std::vector<int>>);
    }

    {
        constexpr auto parser = *cu >> string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::string, std::string>>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      tuple<std::string, std::string>>);
    }

    {
        constexpr auto parser = char_ >> string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }

    {
        constexpr auto parser = char_ >> ' ' >> string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }

    {
        constexpr auto parser = char_ >> ' ' >> string("str") >> ' ' >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }

    {
        constexpr auto parser =
            eps >> char_ >> ' ' >> string("str") >> ' ' >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
        static_assert(std::is_same_v<
                      attribute_t<
                          decltype(BOOST_PARSER_SUBRANGE(first, last)),
                          decltype(parser)>,
                      std::string>);
    }
}
