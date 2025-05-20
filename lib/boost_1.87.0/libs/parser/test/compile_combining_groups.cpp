// Copyright (C) 2024 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>


using namespace boost::parser;

void compile_combining_groups()
{
    // p >> seq_parser
    {
        constexpr auto parser = char_ >> (char_ >> char_);
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<0>, llong<0>, llong<0>>>);
    }
    {
        constexpr auto parser = char_ >> merge[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<0>, llong<1>, llong<1>>>);
    }
    {
        constexpr auto parser = char_ >> separate[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<0>, llong<-1>, llong<-1>>>);
    }



    // seq_parser >> p
    {
        constexpr auto parser = (char_ >> char_) >> char_;
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<0>, llong<0>, llong<0>>>);
    }
    {
        constexpr auto parser = merge[char_ >> char_] >> char_;
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<1>, llong<1>, llong<0>>>);
    }
    {
        constexpr auto parser = separate[char_ >> char_] >> char_;
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<-1>, llong<-1>, llong<0>>>);
    }



    // Apple Clang before 11 produces tuples containing const llong<N>s,
    // causing these tests to spuriously fail.
#if !(defined(__apple_build_version__) && __clang_major__ <= 11)
    // seq_parser >> seq_parser
    {
        constexpr auto parser = (char_ >> char_) >> (char_ >> char_);
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<0>, llong<0>, llong<0>, llong<0>>>);
    }
    {
        constexpr auto parser = merge[char_ >> char_] >> (char_ >> char_);
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<1>, llong<1>, llong<0>, llong<0>>>);
    }
    {
        constexpr auto parser = separate[char_ >> char_] >> (char_ >> char_);
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<-1>, llong<-1>, llong<0>, llong<0>>>);
    }
    {
        constexpr auto parser = (char_ >> char_) >> merge[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<0>, llong<0>, llong<1>, llong<1>>>);
    }
    {
        constexpr auto parser = merge[char_ >> char_] >> merge[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<1>, llong<1>, llong<2>, llong<2>>>);
    }
    {
        constexpr auto parser = separate[char_ >> char_] >> merge[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<-1>, llong<-1>, llong<1>, llong<1>>>);
    }
    {
        constexpr auto parser = (char_ >> char_) >> separate[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<0>, llong<0>, llong<-1>, llong<-1>>>);
    }
    {
        constexpr auto parser = merge[char_ >> char_] >> separate[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<1>, llong<1>, llong<-1>, llong<-1>>>);
    }
    {
        constexpr auto parser = separate[char_ >> char_] >> separate[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<llong<-1>, llong<-1>, llong<-1>, llong<-1>>>);
    }



    // more compilcated cases
    {
        constexpr auto parser = merge[char_ >> char_] >>
                                merge[char_ >> char_] >>
                                merge[char_ >> char_] >> merge[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<
                          llong<1>,
                          llong<1>,
                          llong<2>,
                          llong<2>,
                          llong<3>,
                          llong<3>,
                          llong<4>,
                          llong<4>>>);
    }
    {
        constexpr auto parser = merge[char_ >> char_] >> char_ >> char_ >>
                                merge[char_ >> char_] >> merge[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<
                          llong<1>,
                          llong<1>,
                          llong<0>,
                          llong<0>,
                          llong<2>,
                          llong<2>,
                          llong<3>,
                          llong<3>>>);
    }
    {
        constexpr auto parser = merge[char_ >> char_] >> char_ >> char_ >>
                                separate[char_ >> char_] >> merge[char_ >> char_];
        using combinging_groups_t =
            typename decltype(parser.parser_)::combining_groups;
        static_assert(std::is_same_v<
                      combinging_groups_t,
                      tuple<
                          llong<1>,
                          llong<1>,
                          llong<0>,
                          llong<0>,
                          llong<-1>,
                          llong<-1>,
                          llong<2>,
                          llong<2>>>);
    }
#endif
}
