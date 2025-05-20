/**
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/parser.hpp>
#include <boost/parser/transcode_view.hpp>

#include <boost/core/lightweight_test.hpp>


namespace bp = boost::parser;

bp::symbols<char> const cu_escapes = {{"t", '\t'}, {"r", '\r'}, {"n", '\n'}};
bp::symbols<char32_t> const cp_escapes = {
    {"t", '\t'}, {"r", '\r'}, {"n", '\n'}};

int main()
{

// basic
{
    constexpr auto parser = bp::quoted_string;

    {
        auto result = bp::parse("", parser, bp::ws);
        BOOST_TEST(!result);
    }

    {
        auto result = bp::parse(R"("foo")", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo");
    }

    {
        auto result = bp::parse(R"("foo\\")", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo\\");
    }

    {
        auto result = bp::parse(R"("\"foo\"")", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "\"foo\"");
    }
}

// different_char
{
    constexpr auto parser = bp::quoted_string('\'');

    {
        auto result = bp::parse("", parser, bp::ws);
        BOOST_TEST(!result);
    }

    {
        auto result = bp::parse(R"('foo')", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo");
    }

    {
        auto result = bp::parse(R"('foo\\')", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo\\");
    }

    {
        auto result = bp::parse(R"('\'foo\'')", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "'foo'");
    }
}

// different_char_with_escapes
{
    {
        auto parser = bp::quoted_string('\'', cu_escapes);

        {
            auto result = bp::parse("", parser, bp::ws);
            BOOST_TEST(!result);
        }

        {
            auto result = bp::parse(R"('foo\t')", parser, bp::ws);
            BOOST_TEST(result);
            BOOST_TEST(*result == "foo\t");
        }

        {
            auto result = bp::parse(R"('foo\x')", parser, bp::ws);
            BOOST_TEST(!result);
        }
    }
    {
        auto parser = bp::quoted_string('\'', cp_escapes);

        {
            auto result = bp::parse("", parser, bp::ws);
            BOOST_TEST(!result);
        }

        {
            auto result = bp::parse(R"('\tfoo')", parser, bp::ws);
            BOOST_TEST(result);
            BOOST_TEST(*result == "\tfoo");
        }

        {
            auto result = bp::parse(R"('f\xoo')", parser, bp::ws);
            BOOST_TEST(!result);
        }
    }
}

// char_set
{
    constexpr auto parser = bp::quoted_string("'\"");

    {
        auto result = bp::parse("", parser, bp::ws);
        BOOST_TEST(!result);
    }

    {
        BOOST_TEST(!bp::parse(R"('foo")", parser, bp::ws));
        BOOST_TEST(!bp::parse(R"("foo')", parser, bp::ws));
    }

    {
        auto result = bp::parse(R"('foo')", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo");
    }
    {
        auto result = bp::parse(R"("foo")", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo");
    }

    {
        auto result = bp::parse(R"('foo\\')", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo\\");
    }
    {
        auto result = bp::parse(R"("foo\\")", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "foo\\");
    }

    {
        auto result = bp::parse(R"('\'foo\'')", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "'foo'");
    }
    {
        auto result = bp::parse(R"("\"foo\"")", parser, bp::ws);
        BOOST_TEST(result);
        BOOST_TEST(*result == "\"foo\"");
    }

    {
        // Can't escape arbitrary characters, only backslash and the quote
        // character.
        BOOST_TEST(!bp::parse(R"("\'foo")", parser, bp::ws));
    }
}

// char_set_with_escapes
{
    {
        auto parser = bp::quoted_string("'\"", cu_escapes);

        {
            auto result = bp::parse("", parser, bp::ws);
            BOOST_TEST(!result);
        }

        {
            BOOST_TEST(!bp::parse(R"('foo")", parser, bp::ws));
            BOOST_TEST(!bp::parse(R"("foo')", parser, bp::ws));
        }

        {
            auto result = bp::parse(R"('foo\t')", parser, bp::ws);
            BOOST_TEST(result);
            BOOST_TEST(*result == "foo\t");
        }
        {
            auto result = bp::parse(R"("\tfoo")", parser, bp::ws);
            BOOST_TEST(result);
            BOOST_TEST(*result == "\tfoo");
        }

        {
            auto result = bp::parse(R"('foo\x')", parser, bp::ws);
            BOOST_TEST(!result);
        }
    }
    {
        auto parser = bp::quoted_string("'\"", cp_escapes);

        {
            auto result = bp::parse("", parser, bp::ws);
            BOOST_TEST(!result);
        }

        {
            BOOST_TEST(!bp::parse(R"('foo")", parser, bp::ws));
            BOOST_TEST(!bp::parse(R"("foo')", parser, bp::ws));
        }

        {
            auto result = bp::parse(R"('foo\t')", parser, bp::ws);
            BOOST_TEST(result);
            BOOST_TEST(*result == "foo\t");
        }
        {
            auto result = bp::parse(R"("\tfoo")", parser, bp::ws);
            BOOST_TEST(result);
            BOOST_TEST(*result == "\tfoo");
        }

        {
            auto result = bp::parse(R"('foo\x')", parser, bp::ws);
            BOOST_TEST(!result);
        }
    }
}

// doc_examples
{
    //[ quoted_string_example_1_2
    namespace bp = boost::parser;

    auto result1 = bp::parse("\"some text\"", bp::quoted_string, bp::ws);
    assert(result1);
    std::cout << *result1 << "\n"; // Prints: some text

    auto result2 =
        bp::parse("\"some \\\"text\\\"\"", bp::quoted_string, bp::ws);
    assert(result2);
    std::cout << *result2 << "\n"; // Prints: some "text"
    //]

    //[ quoted_string_example_3
    auto result3 = bp::parse("!some text!", bp::quoted_string('!'), bp::ws);
    assert(result3);
    std::cout << *result3 << "\n"; // Prints: some text
    //]

    //[ quoted_string_example_4
    auto result4 = bp::parse("'some text'", bp::quoted_string("'\""), bp::ws);
    assert(result4);
    std::cout << *result4 << "\n"; // Prints: some text
    //]

    //[ quoted_string_example_5
    // the c++ simple escapes
    bp::symbols<char> const escapes = {
        {"'", '\''},
        {"?", '\?'},
        {"a", '\a'},
        {"b", '\b'},
        {"f", '\f'},
        {"n", '\n'},
        {"r", '\r'},
        {"t", '\t'},
        {"v", '\v'}};
    auto result5 =
        bp::parse("\"some text\r\"", bp::quoted_string('"', escapes), bp::ws);
    assert(result5);
    std::cout << *result5 << "\n"; // Prints (with a CRLF newline): some text
    //]
}

return boost::report_errors();
}
