/**
 *   Copyright (C) 2024 T. Zachary Laine
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

#include <set>

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;


void github_issue_36()
{
    namespace bp = boost::parser;

    auto id = bp::lexeme[+bp::char_(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz0123456789")];
    auto ids = +id;

    std::string str;
    std::vector<std::string> vec;

    bp::parse("id1", id, bp::ws, str); // (1)
    BOOST_TEST(str == "id1");
    str.clear();
    bp::parse("id1 id2", ids, bp::ws, vec); // (2)
    BOOST_TEST(vec == std::vector<std::string>({"id1", "id2"}));

    // Intentionally ill-formed.
    // bp::parse("i1 i2", ids, bp::ws, str);   // (3)
}

namespace issue_50 {
    struct X
    {
        char a;
        int b;

        bool operator<(X rhs) const { return a < rhs.a; }
    };

    struct Y
    {
        std::vector<X> x;
        int c;
    };

    struct Y2
    {
        std::set<X> x;
        int c;
    };

    static_assert(
        boost::parser::detail::is_struct_compatible<
            Y,
            boost::parser::
                tuple<std::vector<boost::parser::tuple<char, int>>, int>>());

    static_assert(
        boost::parser::detail::is_struct_compatible<
            Y2,
            boost::parser::
                tuple<std::vector<boost::parser::tuple<char, int>>, int>>());
}

void github_issue_50()
{
    using namespace issue_50;

    namespace bp = boost::parser;

    {
        auto parse_x = bp::char_ >> bp::int_;
        auto parse_y = +parse_x >> bp::int_;

        Y y;
        auto b = bp::parse("d 3 4", parse_y, bp::ws, y);
        BOOST_TEST(b);

        BOOST_TEST(y.x[0].a == 'd');
        BOOST_TEST(y.x[0].b == 3);
        BOOST_TEST(y.c == 4);
    }

    {
        auto parse_x = bp::char_ >> bp::int_;
        auto parse_y = +parse_x >> bp::int_;

        Y2 y;
        auto b = bp::parse("d 3 4", parse_y, bp::ws, y);
        BOOST_TEST(b);

        BOOST_TEST(y.x.begin()->a == 'd');
        BOOST_TEST(y.x.begin()->b == 3);
        BOOST_TEST(y.c == 4);
    }
}

namespace issue_52 {
    struct X
    {
        char a;
        int b;
    };

    struct Y
    {
        std::vector<X> x;
        int c;
    };

    struct Z
    {
        std::vector<Y> y;
        int d;
    };

    struct W
    {
        std::vector<Z> z;
        int e;
    };
}

void github_issue_52()
{
    using namespace issue_52;

    namespace bp = boost::parser;
    auto parse_x = bp::char_ >> bp::int_;
    auto parse_y = +parse_x >> bp::int_;
    auto parse_z = +parse_y >> bp::char_;
    auto parse_w = +parse_z >> bp::int_;

    {
        Z z;
        auto b = bp::parse("d 2 3 c", parse_z, bp::ws, z);
        BOOST_TEST(b);

        BOOST_TEST(z.y[0].x[0].a == 'd');
        BOOST_TEST(z.y[0].x[0].b == 2);
        BOOST_TEST(z.y[0].c == 3);
        BOOST_TEST(z.d == 'c');
    }
    {
        W w;
        auto b = bp::parse("d 2 3 c 4", parse_w, bp::ws, w);
        BOOST_TEST(b);

        BOOST_TEST(w.z[0].y[0].x[0].a == 'd');
        BOOST_TEST(w.z[0].y[0].x[0].b == 2);
        BOOST_TEST(w.z[0].y[0].c == 3);
        BOOST_TEST(w.z[0].d == 'c');
        BOOST_TEST(w.e == 4);
    }
}

void github_issue_78()
{
    namespace bp = boost::parser;
    std::vector<int> result;
    auto b = bp::parse("3 4 c", +bp::int_, bp::ws, result);
    BOOST_TEST(!b);
    BOOST_TEST(result.empty());
}

namespace issue_90 { namespace parser {
    const auto string =
        ('"' >> boost::parser::lexeme[*(boost::parser::char_ - '"')]) > '"';
    const auto specifier = string > ':' > string;
}}

void github_issue_90()
{
    using namespace issue_90;
    namespace bp = boost::parser;

    std::string input = R"( "dd" : "2" )";
    std::pair<std::string, std::string> result;

    auto b = bp::parse(input, parser::specifier, bp::ws, result);
    BOOST_TEST(b);
    BOOST_TEST(result.first == "dd");
    BOOST_TEST(result.second == "2");
}

namespace github_issue_125_ {
    namespace bp = boost::parser;
    constexpr bp::
        rule<struct replacement_field_rule, std::optional<unsigned short>>
            replacement_field = "replacement_field";
    constexpr auto replacement_field_def = bp::lit('{') >> -bp::ushort_;
    BOOST_PARSER_DEFINE_RULES(replacement_field);
}

void github_issue_125()
{
    namespace bp = boost::parser;
    using namespace github_issue_125_;

    unsigned short integer_found = 99;
    auto print_repl_field = [&](auto & ctx) {
        const std::optional<unsigned short> & val = bp::_attr(ctx);
        if (val)
            integer_found = *val;
        else
            integer_found = 77;
    };

    {
        integer_found = 99;
        auto result = bp::parse("{9", replacement_field[print_repl_field]);
        BOOST_TEST(result);
        BOOST_TEST(integer_found == 9);
    }
    {
        integer_found = 99;
        auto result = bp::parse("{", replacement_field[print_repl_field]);
        BOOST_TEST(result);
        BOOST_TEST(integer_found == 77);
    }
}


int main()
{
    github_issue_36();
    github_issue_50();
    github_issue_52();
    github_issue_78();
    github_issue_90();
    github_issue_125();
    return boost::report_errors();
}
