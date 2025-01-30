// Copyright (C) 2018 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>


using namespace boost::parser;

void compile_seq_attribute()
{
    char const chars[] = "";
    auto first = std::begin(chars);
    auto const last = std::end(chars);

    // scalar and eps
    {
        constexpr auto parser = int_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<int>>);
    }
    {
        constexpr auto parser = eps >> int_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<int>>);
    }

    // scalar >> scalar
    {
        constexpr auto parser = char_ >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = eps >> char_ >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = char_ >> eps >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = char_ >> char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = int_ >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<tuple<int, char>>>);
    }
    {
        constexpr auto parser = eps >> int_ >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<tuple<int, char>>>);
    }
    {
        constexpr auto parser = int_ >> eps >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<tuple<int, char>>>);
    }
    {
        constexpr auto parser = int_ >> char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<tuple<int, char>>>);
    }

    // -scalar >> -scalar
    {
        constexpr auto parser = -char_ >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<
                          tuple<std::optional<char>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = eps >> -char_ >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<
                          tuple<std::optional<char>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = -char_ >> eps >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<
                          tuple<std::optional<char>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = -char_ >> -char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<
                          tuple<std::optional<char>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = -int_ >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<tuple<std::optional<int>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = eps >> -int_ >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<tuple<std::optional<int>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = -int_ >> eps >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<tuple<std::optional<int>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = -int_ >> -char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<tuple<std::optional<int>, std::optional<char>>>>);
    }

    // seq<T> >> seq<T>
    {
        constexpr auto parser = *char_ >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::string, std::string>>>);
    }
    {
        constexpr auto parser = eps >> *char_ >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::string, std::string>>>);
    }
    {
        constexpr auto parser = *char_ >> eps >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::string, std::string>>>);
    }
    {
        constexpr auto parser = *char_ >> *char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::string, std::string>>>);
    }
    {
        constexpr auto parser = *string("str") >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<
                          std::vector<std::string>,
                          std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = eps >> *string("str") >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<
                          std::vector<std::string>,
                          std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = *string("str") >> eps >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<
                          std::vector<std::string>,
                          std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = *string("str") >> *string("str") >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<
                          std::vector<std::string>,
                          std::vector<std::string>>>>);
    }

    // seq<T> >> seq<U>
    {
        constexpr auto parser = *char_ >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::vector<std::string>>>);
    }
    {
        constexpr auto parser = eps >> *char_ >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::vector<std::string>>>);
    }
    {
        constexpr auto parser = *(char_ - "str") >> eps >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::vector<std::string>>>);
    }
    {
        constexpr auto parser = *(char_ >> eps) >> eps >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::vector<std::string>>>);
    }
    {
        constexpr auto parser = *char_ >> eps >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::vector<std::string>>>);
    }
    {
        constexpr auto parser = *char_ >> *string("str") >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<attr_t, std::optional<std::vector<std::string>>>);
    }

    // seq<T> >> T
    {
        constexpr auto parser = *char_ >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = eps >> *char_ >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = *char_ >> eps >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = *char_ >> char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = *string("str") >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::vector<std::string>, char>>>);
    }
    {
        constexpr auto parser = eps >> *string("str") >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::vector<std::string>, char>>>);
    }
    {
        constexpr auto parser = *string("str") >> eps >> char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::vector<std::string>, char>>>);
    }
    {
        constexpr auto parser = *string("str") >> char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<std::vector<std::string>, char>>>);
    }

    // T >> seq<T>
    {
        constexpr auto parser = char_ >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = eps >> char_ >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = char_ >> eps >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = char_ >> *char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = char_ >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<char, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = eps >> char_ >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<char, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = char_ >> eps >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<char, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = char_ >> *string("str") >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<char, std::vector<std::string>>>>);
    }

    // seq<T> >> std::optional<T>
    {
        constexpr auto parser = *char_ >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = eps >> *char_ >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = *char_ >> eps >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = *char_ >> -char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = *string("str") >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::vector<std::string>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = eps >> *string("str") >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::vector<std::string>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = *string("str") >> eps >> -char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::vector<std::string>, std::optional<char>>>>);
    }
    {
        constexpr auto parser = *string("str") >> -char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::vector<std::string>, std::optional<char>>>>);
    }

    // std::optional<T> >> seq<T>
    {
        constexpr auto parser = -char_ >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = eps >> -char_ >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = -char_ >> eps >> *char_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = -char_ >> *char_ >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<attr_t, std::optional<std::string>>);
    }
    {
        constexpr auto parser = -char_ >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = eps >> -char_ >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = -char_ >> eps >> *string("str");
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = -char_ >> *string("str") >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }

    // seq grouping
    {
        constexpr auto parser = (-char_ >> *string("str")) >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = -char_ >> (*string("str") >> eps);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = eps >> (-char_ >> *string("str") >> eps);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = eps >> (-char_ >> *string("str") >> eps) >> eps;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(
            std::is_same_v<
                attr_t,
                std::optional<
                    tuple<std::optional<char>, std::vector<std::string>>>>);
    }
    {
        constexpr auto parser = int_ >> string("str") >> double_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<int, std::string, double>>>);
    }
    {
        constexpr auto parser = (int_ >> string("str")) >> double_;
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<int, std::string, double>>>);
    }
    {
        constexpr auto parser = int_ >> (string("str") >> double_);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<int, std::string, double>>>);
    }
    {
        constexpr auto parser = (int_ >> string("str")) >> (double_ >> int_);
        using attr_t = decltype(prefix_parse(first, last, parser));
        static_assert(std::is_same_v<
                      attr_t,
                      std::optional<tuple<int, std::string, double, int>>>);
    }
}
