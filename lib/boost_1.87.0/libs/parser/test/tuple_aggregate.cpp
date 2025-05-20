/**
 *   Copyright (C) 2023 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/parser.hpp>

#if __has_include(<boost/optional/optional.hpp>)
#define TEST_BOOST_OPTIONAL 1
#include <boost/optional/optional.hpp>
#else
#define TEST_BOOST_OPTIONAL 0
#endif

#if __has_include(<boost/variant/variant.hpp>)
#define TEST_BOOST_VARIANT 1
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#else
#define TEST_BOOST_VARIANT 0
#endif

#if __has_include(<boost/variant2/variant.hpp>)
#define TEST_BOOST_VARIANT2 1
#include <boost/variant2/variant.hpp>
#else
#define TEST_BOOST_VARIANT2 0
#endif

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;

#if TEST_BOOST_OPTIONAL
template<typename T>
constexpr bool boost::parser::enable_optional<boost::optional<T>> = true;
#endif


struct s0_rule_tag
{};
struct s1_rule_a_tag
{};
struct s1_rule_b_tag
{};
struct s2_rule_a_tag
{};
struct s2_rule_b_tag
{};


struct s0
{
    int i_;
    std::string str_;
    std::vector<int> vec_;
};

using s0_tuple = tuple<int, std::string, std::vector<int>>;

auto s0_parser = "s0" >> int_ >> lexeme[+(char_ - ' ')] >> *int_;

callback_rule<s0_rule_tag, s0> s0_rule = "s0_rule";
auto s0_rule_def = s0_parser;

struct s1
{
    int i_;
    std::variant<std::string, double> str_;
    std::vector<int> vec_;
};
#if TEST_BOOST_VARIANT
struct s1_boost_variant
{
    int i_;
    boost::variant<std::string, double> str_;
    std::vector<int> vec_;
};
#endif
#if TEST_BOOST_VARIANT2
struct s1_boost_variant2
{
    int i_;
    boost::variant2::variant<std::string, double> str_;
    std::vector<int> vec_;
};
#endif

using s1_tuple =
    tuple<int, std::variant<std::string, double>, std::vector<int>>;

auto s1_parser_a = "s1" >> int_ >> lexeme[+(char_ - ' ')] >> *int_;
auto s1_parser_b = "s1" >> int_ >> double_ >> *int_;

callback_rule<s1_rule_a_tag, s1> s1_rule_a = "s1_rule_a";
callback_rule<s1_rule_b_tag, s1> s1_rule_b = "s1_rule_b";
auto s1_rule_a_def = s1_parser_a;
auto s1_rule_b_def = s1_parser_b;
#if TEST_BOOST_VARIANT
callback_rule<struct s1_bv_rule_a_tag, s1_boost_variant>
    s1_boost_variant_rule_a = "s1_rule_a";
callback_rule<struct s1_bv_rule_b_tag, s1_boost_variant>
    s1_boost_variant_rule_b = "s1_rule_b";
auto s1_boost_variant_rule_a_def = s1_parser_a;
auto s1_boost_variant_rule_b_def = s1_parser_b;
BOOST_PARSER_DEFINE_RULES(s1_boost_variant_rule_a, s1_boost_variant_rule_b);
#endif
#if TEST_BOOST_VARIANT2
callback_rule<struct s1_bv2_rule_a_tag, s1_boost_variant2>
    s1_boost_variant2_rule_a = "s1_rule_a";
callback_rule<struct s1_bv2_rule_b_tag, s1_boost_variant2>
    s1_boost_variant2_rule_b = "s1_rule_b";
auto s1_boost_variant2_rule_a_def = s1_parser_a;
auto s1_boost_variant2_rule_b_def = s1_parser_b;
BOOST_PARSER_DEFINE_RULES(s1_boost_variant2_rule_a, s1_boost_variant2_rule_b);
#endif

struct s2
{
    int i_;
    std::string str_;
    std::optional<std::vector<int>> vec_;
};
#if TEST_BOOST_OPTIONAL
struct s2_boost_optional
{
    int i_;
    std::string str_;
    boost::optional<std::vector<int>> vec_;
};
#endif

using s2_tuple = tuple<int, std::string, std::optional<std::vector<int>>>;

auto s2_parser_a = "s2" >> int_ >> lexeme[+(char_ - ' ')] >> *int_;
auto s2_parser_b = "s2" >> int_ >> lexeme[+(char_ - ' ')] >> -+int_;

callback_rule<s2_rule_a_tag, s2> s2_rule_a = "s2_rule_a";
callback_rule<s2_rule_b_tag, s2> s2_rule_b = "s2_rule_b";
auto s2_rule_a_def = s2_parser_a;
auto s2_rule_b_def = s2_parser_b;
#if TEST_BOOST_OPTIONAL
callback_rule<struct s2_bo_rule_a_tag, s2_boost_optional>
    s2_boost_optional_rule_a = "s2_rule_a";
callback_rule<struct s2_bo_rule_b_tag, s2_boost_optional>
    s2_boost_optional_rule_b = "s2_rule_b";
auto s2_boost_optional_rule_a_def = s2_parser_a;
auto s2_boost_optional_rule_b_def = s2_parser_b;
BOOST_PARSER_DEFINE_RULES(s2_boost_optional_rule_a, s2_boost_optional_rule_b);
#endif

BOOST_PARSER_DEFINE_RULES(s0_rule, s1_rule_a, s1_rule_b, s2_rule_a, s2_rule_b);

struct s0_like
{
    int64_t i_;
    std::string str_;
    std::vector<int> vec_;
};

struct callbacks_t
{
    void operator()(s0_rule_tag, s0 s) const { s0s.push_back(std::move(s)); }
    void operator()(s1_rule_a_tag, s1 s) const { s1s.push_back(std::move(s)); }
    void operator()(s1_rule_b_tag, s1 s) const { s1s.push_back(std::move(s)); }
    void operator()(s2_rule_a_tag, s2 s) const { s2s.push_back(std::move(s)); }
    void operator()(s2_rule_b_tag, s2 s) const { s2s.push_back(std::move(s)); }

    mutable std::vector<s0> s0s;
    mutable std::vector<s1> s1s;
    mutable std::vector<s2> s2s;
};

int main()
{

// seq_parser_struct_rule
{
    ////////////////////////////////////////////////////////////////////////////
    // Parse-generated attribute.

    {
        std::optional<s0> result = parse("s0 42 text 1 2 3", s0_rule, ws);
        BOOST_TEST(result);
        s0 & struct_ = *result;
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<s1> const result =
            parse("s1 42 text 1 2 3", s1_rule_a, ws);
        BOOST_TEST(result);
        s1 const & struct_ = *result;
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<s1> const result =
            parse("s1 42 13.0 1 2 3", s1_rule_b, ws);
        BOOST_TEST(result);
        s1 const & struct_ = *result;
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<s2> const result =
            parse("s2 42 text 1 2 3", s2_rule_a, ws);
        BOOST_TEST(result);
        s2 const & struct_ = *result;
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<s2> const result =
            parse("s2 42 text 1 2 3", s2_rule_b, ws);
        BOOST_TEST(result);
        s2 const & struct_ = *result;
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }

    // Use the rule as part of a larger parse.
    {
        std::optional<tuple<int, s0>> const result =
            parse("99 s0 42 text 1 2 3", int_ >> s0_rule, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(*result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<tuple<int, s1>> const result =
            parse("99 s1 42 text 1 2 3", int_ >> s1_rule_a, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(*result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<tuple<int, s1>> const result =
            parse("99 s1 42 13.0 1 2 3", int_ >> s1_rule_b, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(*result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<tuple<int, s2>> const result =
            parse("99 s2 42 text 1 2 3", int_ >> s2_rule_a, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(*result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        std::optional<tuple<int, s2>> const result =
            parse("99 s2 42 text 1 2 3", int_ >> s2_rule_b, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(*result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }

    ////////////////////////////////////////////////////////////////////////////
    // Pass attribute to parse.

    {
        s0 struct_;

        BOOST_TEST(parse("s0 42 text 1 2 3", s0_rule, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1 struct_;

        BOOST_TEST(parse("s1 42 text 1 2 3", s1_rule_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1 struct_;

        BOOST_TEST(parse("s1 42 13.0 1 2 3", s1_rule_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#if TEST_BOOST_VARIANT
    {
        s1_boost_variant struct_;

        BOOST_TEST(
            parse("s1 42 text 1 2 3", s1_boost_variant_rule_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).which() == 0u);
        BOOST_TEST(boost::get<std::string>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_boost_variant struct_;

        BOOST_TEST(
            parse("s1 42 13.0 1 2 3", s1_boost_variant_rule_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).which() == 1u);
        BOOST_TEST(boost::get<double>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
#if TEST_BOOST_VARIANT2
    {
        s1_boost_variant2 struct_;

        BOOST_TEST(
            parse("s1 42 text 1 2 3", s1_boost_variant2_rule_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(boost::variant2::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_boost_variant2 struct_;

        BOOST_TEST(
            parse("s1 42 13.0 1 2 3", s1_boost_variant2_rule_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(boost::variant2::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
    {
        s2 struct_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_rule_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2 struct_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_rule_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#if TEST_BOOST_OPTIONAL
    {
        s2_boost_optional struct_;

        BOOST_TEST(
            parse("s2 42 text 1 2 3", s2_boost_optional_rule_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2_boost_optional struct_;

        BOOST_TEST(
            parse("s2 42 text 1 2 3", s2_boost_optional_rule_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif

    // Use the rule as part of a larger parse.
    {
        tuple<int, s0> result;

        BOOST_TEST(parse("99 s0 42 text 1 2 3", int_ >> s0_rule, ws, result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        tuple<int, s1> result;

        BOOST_TEST(parse("99 s1 42 text 1 2 3", int_ >> s1_rule_a, ws, result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        tuple<int, s1> result;

        BOOST_TEST(parse("99 s1 42 13.0 1 2 3", int_ >> s1_rule_b, ws, result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#if TEST_BOOST_VARIANT
    {
        tuple<int, s1_boost_variant> result;

        BOOST_TEST(parse(
            "99 s1 42 text 1 2 3",
            int_ >> s1_boost_variant_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).which() == 0u);
        BOOST_TEST(boost::get<std::string>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        tuple<int, s1_boost_variant> result;

        BOOST_TEST(parse(
            "99 s1 42 13.0 1 2 3",
            int_ >> s1_boost_variant_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).which() == 1u);
        BOOST_TEST(boost::get<double>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
#if TEST_BOOST_VARIANT2
    {
        tuple<int, s1_boost_variant2> result;

        BOOST_TEST(parse(
            "99 s1 42 text 1 2 3",
            int_ >> s1_boost_variant2_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(boost::variant2::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        tuple<int, s1_boost_variant2> result;

        BOOST_TEST(parse(
            "99 s1 42 13.0 1 2 3",
            int_ >> s1_boost_variant2_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(boost::variant2::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
    {
        tuple<int, s2> result;

        BOOST_TEST(parse("99 s2 42 text 1 2 3", int_ >> s2_rule_a, ws, result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        tuple<int, s2> result;

        BOOST_TEST(parse("99 s2 42 text 1 2 3", int_ >> s2_rule_b, ws, result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#if TEST_BOOST_OPTIONAL
    {
        tuple<int, s2_boost_optional> result;

        BOOST_TEST(parse(
            "99 s2 42 text 1 2 3",
            int_ >> s2_boost_optional_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        tuple<int, s2_boost_optional> result;

        BOOST_TEST(parse(
            "99 s2 42 text 1 2 3",
            int_ >> s2_boost_optional_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto struct_ = get(result, llong<1>{});
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
}

// repeated_seq_parser_struct_rule
{
    ////////////////////////////////////////////////////////////////////////////
    // Parse-generated attribute.

    {
        std::optional<std::vector<s0>> result =
            parse("s0 42 text 1 2 3 s0 41 texty 1 3 2", *s0_rule, ws);
        BOOST_TEST(result);
        std::vector<s0> & structs_ = *result;
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<std::vector<s1>> const result =
            parse("s1 42 text 1 2 3 s1 41 texty 1 3 2", *s1_rule_a, ws);
        BOOST_TEST(result);
        std::vector<s1> const & structs_ = *result;
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<std::vector<s1>> const result =
            parse("s1 42 13.0 1 2 3 s1 41 12.0 1 3 2", *s1_rule_b, ws);
        BOOST_TEST(result);
        std::vector<s1> const & structs_ = *result;
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<std::vector<s2>> const result =
            parse("s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_rule_a, ws);
        BOOST_TEST(result);
        std::vector<s2> const & structs_ = *result;
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<std::vector<s2>> const result =
            parse("s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_rule_b, ws);
        BOOST_TEST(result);
        std::vector<s2> const & structs_ = *result;
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }

    // Use the rule as part of a larger parse.
    {
        std::optional<tuple<int, std::vector<s0>>> const result =
            parse("99 s0 42 text 1 2 3 s0 41 texty 1 3 2", int_ >> *s0_rule, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(*result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<tuple<int, std::vector<s1>>> const result =
            parse("99 s1 42 text 1 2 3 s1 41 texty 1 3 2", int_ >> *s1_rule_a, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(*result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<tuple<int, std::vector<s1>>> const result =
            parse("99 s1 42 13.0 1 2 3 s1 41 12.0 1 3 2", int_ >> *s1_rule_b, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(*result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<tuple<int, std::vector<s2>>> const result =
            parse("99 s2 42 text 1 2 3 s2 41 texty 1 3 2", int_ >> *s2_rule_a, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(*result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::optional<tuple<int, std::vector<s2>>> const result =
            parse("99 s2 42 text 1 2 3 s2 41 texty 1 3 2", int_ >> *s2_rule_b, ws);
        auto i_ = get(*result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(*result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }

    ////////////////////////////////////////////////////////////////////////////
    // Pass attribute to parse.

    {
        std::vector<s0> structs_;

        BOOST_TEST(parse(
            "s0 42 text 1 2 3 s0 41 texty 1 3 2", *s0_rule, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1> structs_;

        BOOST_TEST(parse(
            "s1 42 text 1 2 3 s1 41 texty 1 3 2", *s1_rule_a, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1> structs_;

        BOOST_TEST(parse(
            "s1 42 13.0 1 2 3 s1 41 12.0 1 3 2", *s1_rule_b, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#if TEST_BOOST_VARIANT
    {
        std::vector<s1_boost_variant> structs_;

        BOOST_TEST(parse(
            "s1 42 text 1 2 3 s1 41 texty 1 3 2",
            *s1_boost_variant_rule_a,
            ws,
            structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).which() == 0u);
        BOOST_TEST(
            boost::get<std::string>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).which() == 0u);
        BOOST_TEST(
            boost::get<std::string>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1_boost_variant> structs_;

        BOOST_TEST(parse(
            "s1 42 13.0 1 2 3 s1 41 12.0 1 3 2",
            *s1_boost_variant_rule_b,
            ws,
            structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).which() == 1u);
        BOOST_TEST(boost::get<double>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).which() == 1u);
        BOOST_TEST(boost::get<double>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#endif
#if TEST_BOOST_VARIANT2
    {
        std::vector<s1_boost_variant2> structs_;

        BOOST_TEST(parse(
            "s1 42 text 1 2 3 s1 41 texty 1 3 2",
            *s1_boost_variant2_rule_a,
            ws,
            structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(
            boost::variant2::get<0>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(
            boost::variant2::get<0>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1_boost_variant2> structs_;

        BOOST_TEST(parse(
            "s1 42 13.0 1 2 3 s1 41 12.0 1 3 2",
            *s1_boost_variant2_rule_b,
            ws,
            structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(boost::variant2::get<1>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(boost::variant2::get<1>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#endif
    {
        std::vector<s2> structs_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_rule_a, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2> structs_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_rule_b, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#if TEST_BOOST_OPTIONAL
    {
        std::vector<s2_boost_optional> structs_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2",
            *s2_boost_optional_rule_a,
            ws,
            structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2_boost_optional> structs_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2",
            *s2_boost_optional_rule_b,
            ws,
            structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#endif

    // Use the rule as part of a larger parse.
    {
        tuple<int, std::vector<s0>> result;

        BOOST_TEST(parse(
            "99 s0 42 text 1 2 3 s0 41 texty 1 3 2",
            int_ >> *s0_rule,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        tuple<int, std::vector<s1>> result;

        BOOST_TEST(parse(
            "99 s1 42 text 1 2 3 s1 41 texty 1 3 2",
            int_ >> *s1_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        tuple<int, std::vector<s1>> result;

        BOOST_TEST(parse(
            "99 s1 42 13.0 1 2 3 s1 41 12.0 1 3 2",
            int_ >> *s1_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#if TEST_BOOST_VARIANT
    {
        tuple<int, std::vector<s1_boost_variant>> result;

        BOOST_TEST(parse(
            "99 s1 42 text 1 2 3 s1 41 texty 1 3 2",
            int_ >> *s1_boost_variant_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).which() == 0u);
        BOOST_TEST(
            boost::get<std::string>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).which() == 0u);
        BOOST_TEST(
            boost::get<std::string>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        tuple<int, std::vector<s1_boost_variant>> result;

        BOOST_TEST(parse(
            "99 s1 42 13.0 1 2 3 s1 41 12.0 1 3 2",
            int_ >> *s1_boost_variant_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).which() == 1u);
        BOOST_TEST(boost::get<double>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).which() == 1u);
        BOOST_TEST(boost::get<double>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#endif
#if TEST_BOOST_VARIANT2
    {
        tuple<int, std::vector<s1_boost_variant2>> result;

        BOOST_TEST(parse(
            "99 s1 42 text 1 2 3 s1 41 texty 1 3 2",
            int_ >> *s1_boost_variant2_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(
            boost::variant2::get<0>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(
            boost::variant2::get<0>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        tuple<int, std::vector<s1_boost_variant2>> result;

        BOOST_TEST(parse(
            "99 s1 42 13.0 1 2 3 s1 41 12.0 1 3 2",
            int_ >> *s1_boost_variant2_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(boost::variant2::get<1>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(boost::variant2::get<1>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#endif
    {
        tuple<int, std::vector<s2>> result;

        BOOST_TEST(parse(
            "99 s2 42 text 1 2 3 s2 41 texty 1 3 2",
            int_ >> *s2_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        tuple<int, std::vector<s2>> result;

        BOOST_TEST(parse(
            "99 s2 42 text 1 2 3 s2 41 texty 1 3 2",
            int_ >> *s2_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#if TEST_BOOST_OPTIONAL
    {
        tuple<int, std::vector<s2_boost_optional>> result;

        BOOST_TEST(parse(
            "99 s2 42 text 1 2 3 s2 41 texty 1 3 2",
            int_ >> *s2_boost_optional_rule_a,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        tuple<int, std::vector<s2_boost_optional>> result;

        BOOST_TEST(parse(
            "99 s2 42 text 1 2 3 s2 41 texty 1 3 2",
            int_ >> *s2_boost_optional_rule_b,
            ws,
            result));
        auto i_ = get(result, llong<0>{});
        BOOST_TEST(i_ == 99);
        auto structs_ = get(result, llong<1>{});
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
#endif
}

// seq_parser_struct_cb_rule
{
    {
        callbacks_t callbacks;
        BOOST_TEST(callback_parse("s0 42 text 1 2 3", s0_rule, ws, callbacks));
        BOOST_TEST(callbacks.s0s.size() == 1u);
        s0 const & struct_ = callbacks.s0s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
            BOOST_TEST(callback_parse("s1 42 text 1 2 3", s1_rule_a, ws, callbacks));
        BOOST_TEST(callbacks.s1s.size() == 1u);
        s1 const & struct_ = callbacks.s1s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
            BOOST_TEST(callback_parse("s1 42 13.0 1 2 3", s1_rule_b, ws, callbacks));
        BOOST_TEST(callbacks.s1s.size() == 1u);
        s1 const & struct_ = callbacks.s1s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
            BOOST_TEST(callback_parse("s2 42 text 1 2 3", s2_rule_a, ws, callbacks));
        BOOST_TEST(callbacks.s2s.size() == 1u);
        s2 const & struct_ = callbacks.s2s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
            BOOST_TEST(callback_parse("s2 42 text 1 2 3", s2_rule_b, ws, callbacks));
        BOOST_TEST(callbacks.s2s.size() == 1u);
        s2 const & struct_ = callbacks.s2s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }

    // Use the rule as part of a larger parse.
    {
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(
            "99 s0 42 text 1 2 3", int_ >> s0_rule, ws, callbacks));
        BOOST_TEST(callbacks.s0s.size() == 1u);
        s0 const & struct_ = callbacks.s0s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(
            "99 s1 42 text 1 2 3", int_ >> s1_rule_a, ws, callbacks));
        BOOST_TEST(callbacks.s1s.size() == 1u);
        s1 const & struct_ = callbacks.s1s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(
            "99 s1 42 13.0 1 2 3", int_ >> s1_rule_b, ws, callbacks));
        BOOST_TEST(callbacks.s1s.size() == 1u);
        s1 const & struct_ = callbacks.s1s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(
            "99 s2 42 text 1 2 3", int_ >> s2_rule_a, ws, callbacks));
        BOOST_TEST(callbacks.s2s.size() == 1u);
        s2 const & struct_ = callbacks.s2s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(
            "99 s2 42 text 1 2 3", int_ >> s2_rule_b, ws, callbacks));
        BOOST_TEST(callbacks.s2s.size() == 1u);
        s2 const & struct_ = callbacks.s2s[0];
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
}

// parse_into_struct
{
    // tuples
    {
        s0_tuple tuple_;

        BOOST_TEST(parse("s0 42 text 1 2 3", s0_parser, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_tuple tuple_;

        BOOST_TEST(parse("s1 42 text 1 2 3", s1_parser_a, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(tuple_, llong<1>{})) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_tuple tuple_;

        BOOST_TEST(parse("s1 42 13.0 1 2 3", s1_parser_b, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(tuple_, llong<1>{})) == 13.0);
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2_tuple tuple_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_parser_a, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2_tuple tuple_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_parser_b, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }

    // structs
    {
        s0 struct_;

        BOOST_TEST(parse("s0 42 text 1 2 3", s0_parser, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s0_like struct_;

        BOOST_TEST(parse("s0 42 text 1 2 3", s0_parser, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1 struct_;

        BOOST_TEST(parse("s1 42 text 1 2 3", s1_parser_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1 struct_;

        BOOST_TEST(parse("s1 42 13.0 1 2 3", s1_parser_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#if TEST_BOOST_VARIANT
    {
        s1_boost_variant struct_;

        BOOST_TEST(parse("s1 42 text 1 2 3", s1_parser_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).which() == 0u);
        BOOST_TEST(boost::get<std::string>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_boost_variant struct_;

        BOOST_TEST(parse("s1 42 13.0 1 2 3", s1_parser_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).which() == 1u);
        BOOST_TEST(boost::get<double>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
#if TEST_BOOST_VARIANT2
    {
        s1_boost_variant2 struct_;

        BOOST_TEST(parse("s1 42 text 1 2 3", s1_parser_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 0u);
        BOOST_TEST(boost::variant2::get<0>(get(struct_, llong<1>{})) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_boost_variant2 struct_;

        BOOST_TEST(parse("s1 42 13.0 1 2 3", s1_parser_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}).index() == 1u);
        BOOST_TEST(boost::variant2::get<1>(get(struct_, llong<1>{})) == 13.0);
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
    {
        s2 struct_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_parser_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2 struct_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_parser_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#if TEST_BOOST_OPTIONAL
    {
        s2_boost_optional struct_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_parser_a, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(*get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2_boost_optional struct_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_parser_b, ws, struct_));
        BOOST_TEST(get(struct_, llong<0>{}) == 42);
        BOOST_TEST(get(struct_, llong<1>{}) == "text");
        BOOST_TEST(*get(struct_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
#endif
}

// repeated_parse_into_struct
{
    // tuples
    {
        std::vector<s0_tuple> tuples_;

        BOOST_TEST(parse(
            "s0 42 text 1 2 3 s0 41 texty 1 3 2", *s0_parser, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1_tuple> tuples_;

        BOOST_TEST(parse(
            "s1 42 text 1 2 3 s1 41 texty 1 3 2", *s1_parser_a, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(tuples_[0], llong<1>{})) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(tuples_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1_tuple> tuples_;

        BOOST_TEST(parse(
            "s1 42 13.0 1 2 3 s1 41 12.0 1 3 2", *s1_parser_b, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(tuples_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(tuples_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2_tuple> tuples_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_parser_a, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2_tuple> tuples_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_parser_b, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }

    // structs
    {
        std::vector<s0> structs_;

        BOOST_TEST(parse(
            "s0 42 text 1 2 3 s0 41 texty 1 3 2", *s0_parser, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s0_like> structs_;

        BOOST_TEST(parse(
            "s0 42 text 1 2 3 s0 41 texty 1 3 2", *s0_parser, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1> structs_;

        BOOST_TEST(parse(
            "s1 42 text 1 2 3 s1 41 texty 1 3 2", *s1_parser_a, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[0], llong<1>{})) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(structs_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1> structs_;

        BOOST_TEST(parse(
            "s1 42 13.0 1 2 3 s1 41 12.0 1 3 2", *s1_parser_b, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(structs_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2> structs_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_parser_a, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2> structs_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_parser_b, ws, structs_));
        BOOST_TEST(get(structs_[0], llong<0>{}) == 42);
        BOOST_TEST(get(structs_[0], llong<1>{}) == "text");
        BOOST_TEST(get(structs_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(structs_[1], llong<0>{}) == 41);
        BOOST_TEST(get(structs_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(structs_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
}

// parse_into_tuple
{
    {
        s0_tuple tuple_;

        BOOST_TEST(parse("s0 42 text 1 2 3", s0_rule, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_tuple tuple_;

        BOOST_TEST(parse("s1 42 text 1 2 3", s1_rule_a, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(tuple_, llong<1>{})) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s1_tuple tuple_;

        BOOST_TEST(parse("s1 42 13.0 1 2 3", s1_rule_b, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(tuple_, llong<1>{})) == 13.0);
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2_tuple tuple_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_rule_a, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
    {
        s2_tuple tuple_;

        BOOST_TEST(parse("s2 42 text 1 2 3", s2_rule_b, ws, tuple_));
        BOOST_TEST(get(tuple_, llong<0>{}) == 42);
        BOOST_TEST(get(tuple_, llong<1>{}) == "text");
        BOOST_TEST(get(tuple_, llong<2>{}) == std::vector<int>({1, 2, 3}));
    }
}

// repeated_parse_into_tuple
{
    {
        std::vector<s0_tuple> tuples_;

        BOOST_TEST(
            parse("s0 42 text 1 2 3 s0 41 texty 1 3 2", *s0_rule, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1_tuple> tuples_;

        BOOST_TEST(parse(
            "s1 42 text 1 2 3 s1 41 texty 1 3 2", *s1_rule_a, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(tuples_[0], llong<1>{})) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}).index() == 0u);
        BOOST_TEST(std::get<0>(get(tuples_[1], llong<1>{})) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s1_tuple> tuples_;

        BOOST_TEST(parse(
            "s1 42 13.0 1 2 3 s1 41 12.0 1 3 2", *s1_rule_b, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(tuples_[0], llong<1>{})) == 13.0);
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}).index() == 1u);
        BOOST_TEST(std::get<1>(get(tuples_[1], llong<1>{})) == 12.0);
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2_tuple> tuples_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_rule_a, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
    {
        std::vector<s2_tuple> tuples_;

        BOOST_TEST(parse(
            "s2 42 text 1 2 3 s2 41 texty 1 3 2", *s2_rule_b, ws, tuples_));
        BOOST_TEST(get(tuples_[0], llong<0>{}) == 42);
        BOOST_TEST(get(tuples_[0], llong<1>{}) == "text");
        BOOST_TEST(get(tuples_[0], llong<2>{}) == std::vector<int>({1, 2, 3}));
        BOOST_TEST(get(tuples_[1], llong<0>{}) == 41);
        BOOST_TEST(get(tuples_[1], llong<1>{}) == "texty");
        BOOST_TEST(get(tuples_[1], llong<2>{}) == std::vector<int>({1, 3, 2}));
    }
}

return boost::report_errors();
}
