/**
 *   Copyright (C) 2020 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/parser.hpp>


// This is unfortunately a visual test.  It is not possible to print most
// parsers without the parser context.  So, just run this and look at the
// output.

using namespace boost::parser;

struct globals_t
{
    int i = 42;
    int i2 = 84;
    unsigned int u = 42u;
    char c = 'g';
    uint32_t cp = U'\u0300';
};

globals_t const globals;

auto i = [](auto & ctx) { return _globals(ctx).i; };
auto i2 = [](auto & ctx) { return _globals(ctx).i2; };
auto u = [](auto & ctx) { return _globals(ctx).u; };
auto c = [](auto & ctx) { return _globals(ctx).c; };
[[maybe_unused]] auto cp_ = [](auto & ctx) { return _globals(ctx).cp; };

auto true_ = [](auto & ctx) { return true; };
auto zero = [](auto & ctx) { return 0; };

auto a = [](auto & ctx) {};

rule<class test_rule> const test_rule = "test_rule";
auto const test_rule_def = char_;
BOOST_PARSER_DEFINE_RULES(test_rule);

#define PARSE(p)                                                               \
    std::cout << #p << ":\n";                                                  \
    parse(str, with_globals(p, globals), trace::on)

#define PARSE_CHAR32(p)                                                        \
    std::cout << #p << ":\n";                                                  \
    parse(u32str, with_globals(p, globals), trace::on)

int main()
{
    std::string const str = "";
    std::string const char_range = "xyz";

    std::cout << "----------------------------------------\n"
              << "| repeat()[]                            |\n"
              << "----------------------------------------\n";

    PARSE(repeat(42)[char_]);
    PARSE(repeat(i)[char_]);
    PARSE(repeat(42, 42)[char_]);
    PARSE(repeat(i, i)[char_]);
    PARSE(repeat(42, 84)[char_]);
    PARSE(repeat(i, i2)[char_]);
    PARSE(repeat(0, Inf)[char_]);
    PARSE(repeat(1, Inf)[char_]);
    PARSE(repeat(2, Inf)[char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator*                             |\n"
              << "----------------------------------------\n";

    PARSE(*char_);
    PARSE(*(*char_ >> char_));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator+                             |\n"
              << "----------------------------------------\n";

    PARSE(+char_);
    PARSE(+(*char_ >> char_));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator%                             |\n"
              << "----------------------------------------\n";

    PARSE(int_ % char_);
    PARSE(int_ % ',');

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator- (unary)                     |\n"
              << "----------------------------------------\n";

    PARSE(-char_);
    PARSE(-(*char_ >> char_));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator|                             |\n"
              << "----------------------------------------\n";

    PARSE(char_ | int_ | float_);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator||                            |\n"
              << "----------------------------------------\n";

    PARSE(char_ || int_ || float_);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator- (binary)                    |\n"
              << "----------------------------------------\n";

    PARSE(char_ - int_);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator>>/operator>                  |\n"
              << "----------------------------------------\n";

#if defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-shift-op-parentheses"
#endif
    PARSE(char_ >> int_ > float_);
#if defined(__clang__)
#pragma GCC diagnostic pop
#endif

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator[]                            |\n"
              << "----------------------------------------\n";

    PARSE(char_[a]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| transform)f_[]                        |\n"
              << "----------------------------------------\n";

    auto f = [](auto x) { return x; };
    PARSE(transform(f)[char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| omit[]                                |\n"
              << "----------------------------------------\n";

    PARSE(omit[char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| lexeme[]                              |\n"
              << "----------------------------------------\n";

    PARSE(lexeme[char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| raw[]                                 |\n"
              << "----------------------------------------\n";

    PARSE(raw[char_]);

#if BOOST_PARSER_USE_CONCEPTS
    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| string_view[]                         |\n"
              << "----------------------------------------\n";

    PARSE(string_view[char_]);
#endif

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| skip[]                                |\n"
              << "----------------------------------------\n";

    PARSE(skip[char_]);
    PARSE(skip(ws)[char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| merge[]                               |\n"
              << "----------------------------------------\n";

    PARSE(merge[char_ >> char_]);
    PARSE(char_ >> merge[char_ >> char_]);
    PARSE(merge[char_ >> char_] >> char_);
    PARSE(merge[char_ >> char_] >> char_ >> char_ >> char_);
    PARSE(merge[char_ >> char_ >> char_ >> char_ >> char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| separate[]                            |\n"
              << "----------------------------------------\n";

    PARSE(separate[char_ >> char_]);
    PARSE(char_ >> separate[char_ >> char_]);
    PARSE(separate[char_ >> char_] >> char_);
    PARSE(separate[char_ >> char_] >> char_ >> char_ >> char_);
    PARSE(separate[char_ >> char_ >> char_ >> char_ >> char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator!                             |\n"
              << "----------------------------------------\n";

    PARSE(!char_);
    PARSE(!(*char_ >> char_));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| operator&                             |\n"
              << "----------------------------------------\n";

    PARSE(&char_);
    PARSE(&(*char_ >> char_));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| symbols<T>                            |\n"
              << "----------------------------------------\n";

    symbols<int> symtab1;
    symbols<float> symtab2 = "integer symbols";
    PARSE(symtab1);
    PARSE(symtab2);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| rules                                 |\n"
              << "----------------------------------------\n";

    PARSE(test_rule);
    PARSE(test_rule.with(0, i));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| eps                                   |\n"
              << "----------------------------------------\n";

    PARSE(eps);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| eoi                                   |\n"
              << "----------------------------------------\n";

    PARSE(eoi);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| attr()                                |\n"
              << "----------------------------------------\n";

    PARSE(attr(42));
    PARSE(attr(i));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| char_                                 |\n"
              << "----------------------------------------\n";

    PARSE(char_);
    PARSE(char_('g'));
    PARSE(char_(c));
    PARSE(char_('g', 'g'));
    PARSE(char_(c, c));
    PARSE(char_("xyz"));
    PARSE(char_(char_range));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| cp                                    |\n"
              << "----------------------------------------\n";

    PARSE(cp);
    PARSE(cp(42u));
    PARSE(cp(u));
    PARSE(cp('g', 'g'));
    PARSE(cp(c, c));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| cu                                    |\n"
              << "----------------------------------------\n";

    PARSE(cu);
    PARSE(cu('g'));
    PARSE(cu(c));
    PARSE(cu('g', 'g'));
    PARSE(cu(c, c));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| lit()                                 |\n"
              << "----------------------------------------\n";

    PARSE(lit('g'));
    PARSE(lit(u8'g'));
    PARSE(lit(U'g'));
    PARSE(lit("g"));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| string()                              |\n"
              << "----------------------------------------\n";

    PARSE(string("h"));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| quoted_string()                       |\n"
              << "----------------------------------------\n";

    PARSE(quoted_string);
    PARSE(quoted_string('\''));
    PARSE(quoted_string("'\""));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| eol                                   |\n"
              << "----------------------------------------\n";

    PARSE(eol);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| ws                                    |\n"
              << "----------------------------------------\n";

    PARSE(ws);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| blank                                 |\n"
              << "----------------------------------------\n";

    PARSE(blank);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| digit                                 |\n"
              << "----------------------------------------\n";

    PARSE(digit);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| control                               |\n"
              << "----------------------------------------\n";

    PARSE(control);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| punct                                 |\n"
              << "----------------------------------------\n";

    PARSE(punct);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| lower                                 |\n"
              << "----------------------------------------\n";

    PARSE(lower);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| upper                                 |\n"
              << "----------------------------------------\n";

    PARSE(upper);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| hex_digit                             |\n"
              << "----------------------------------------\n";

    PARSE(hex_digit);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| bool_                                 |\n"
              << "----------------------------------------\n";

    PARSE(bool_);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| int_                                  |\n"
              << "----------------------------------------\n";

    PARSE(short_);
    PARSE(short_(42));
    PARSE(short_(i));
    PARSE(int_);
    PARSE(int_(42));
    PARSE(int_(i));
    PARSE(long_);
    PARSE(long_(42));
    PARSE(long_(i));
    PARSE(long_long);
    PARSE(long_long(42));
    PARSE(long_long(i));

    parser_interface<int_parser<int, 16, 2, 3>> const hex_int;
    PARSE(hex_int);
    PARSE(hex_int(42));
    PARSE(hex_int(i));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| uint_                                 |\n"
              << "----------------------------------------\n";

    PARSE(bin);
    PARSE(bin(42u));
    PARSE(bin(u));
    PARSE(oct);
    PARSE(oct(42u));
    PARSE(oct(u));
    PARSE(hex);
    PARSE(hex(42u));
    PARSE(hex(u));
    PARSE(ushort_);
    PARSE(ushort_(42u));
    PARSE(ushort_(u));
    PARSE(uint_);
    PARSE(uint_(42u));
    PARSE(uint_(u));
    PARSE(ulong_);
    PARSE(ulong_(42u));
    PARSE(ulong_(u));
    PARSE(ulong_long);
    PARSE(ulong_long(42u));
    PARSE(ulong_long(u));

    parser_interface<uint_parser<unsigned int, 16, 2, 3>> const hex_uint;
    PARSE(hex_uint);
    PARSE(hex_uint(42u));
    PARSE(hex_uint(u));

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| float_                                |\n"
              << "----------------------------------------\n";

    PARSE(float_);
    PARSE(double_);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| if_                                   |\n"
              << "----------------------------------------\n";

    PARSE(if_(true_)[char_]);

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| switch_                               |\n"
              << "----------------------------------------\n";

    PARSE(switch_(0)(0, char_)(1, int_));
    PARSE(switch_(zero)(zero, char_)(i, int_));

    std::u32string const u32str = U"Straße";

    std::cout << "\n\n"
              << "----------------------------------------\n"
              << "| float_ (char32_t)                     |\n"
              << "----------------------------------------\n";

    PARSE_CHAR32(float_);
    PARSE_CHAR32(double_);
}
