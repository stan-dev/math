//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/constant_string_view.hpp>
#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <initializer_list>

#include "format_common.hpp"
#include "test_common/printing.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

namespace boost {
namespace mysql {
namespace test {

// When formatted, outputs the spec to the output context
struct echo_spec
{
};

// Parse returns begin + N
template <std::size_t N>
struct parse_consume
{
};

// Format may call set_error
struct maybe_error
{
    bool do_error;
};

// Only non-const formattable
struct only_nonconst_lvalues
{
    int id;
};

// Formatter has overloads for both const and non-const values
struct const_and_nonconst
{
};

}  // namespace test

template <>
struct formatter<echo_spec>
{
    string_view spec;

    const char* parse(const char* it, const char* end)
    {
        spec = {it, end};
        return end;
    }

    void format(echo_spec, format_context_base& ctx) const { ctx.append_raw(runtime(spec)); }
};

template <std::size_t N>
struct formatter<parse_consume<N>>
{
    const char* parse(const char* it, const char*) { return it + N; }
    void format(parse_consume<N>, format_context_base&) const {}
};

template <>
struct formatter<maybe_error>
{
    const char* parse(const char* it, const char*) { return it; }

    void format(maybe_error v, format_context_base& ctx) const
    {
        if (v.do_error)
            ctx.add_error(client_errc::unformattable_value);
        else
            format_sql_to(ctx, "{}", v.do_error);
    };
};

template <>
struct formatter<only_nonconst_lvalues>
{
    const char* parse(const char* it, const char*) { return it; }
    void format(only_nonconst_lvalues& val, format_context_base& ctx) const
    {
        format_sql_to(ctx, "only_nonconst_lvalues({})", val.id);
    }
};

template <>
struct formatter<const_and_nonconst>
{
    const char* parse(const char* it, const char*) { return it; }
    void format(const_and_nonconst&, format_context_base& ctx) const { ctx.append_raw("nonconst"); }
    void format(const const_and_nonconst&, format_context_base& ctx) const { ctx.append_raw("const"); }
};

}  // namespace mysql
}  // namespace boost

BOOST_AUTO_TEST_SUITE(test_format_sql_custom)

constexpr format_options opts{utf8mb4_charset, true};

// We pass the correct format spec to parse
BOOST_AUTO_TEST_CASE(parse_passed_format_specs)
{
    std::initializer_list<format_arg> named_args{
        {"name", echo_spec{}}
    };

    // No spec
    BOOST_TEST(format_sql(opts, "{}", echo_spec{}) == "");
    BOOST_TEST(format_sql(opts, "{:}", echo_spec{}) == "");
    BOOST_TEST(format_sql(opts, "{0}", echo_spec{}) == "");
    BOOST_TEST(format_sql(opts, "{0:}", echo_spec{}) == "");
    BOOST_TEST(format_sql(opts, "{name}", named_args) == "");
    BOOST_TEST(format_sql(opts, "{name:}", named_args) == "");

    // Single char
    BOOST_TEST(format_sql(opts, "{:k}", echo_spec{}) == "k");
    BOOST_TEST(format_sql(opts, "{0:u}", echo_spec{}) == "u");
    BOOST_TEST(format_sql(opts, "{name:p}", named_args) == "p");

    // Multiple chars
    BOOST_TEST(format_sql(opts, "{:some}", echo_spec{}) == "some");
    BOOST_TEST(format_sql(opts, "{0:chars}", echo_spec{}) == "chars");
    BOOST_TEST(format_sql(opts, "{name:here}", named_args) == "here");

    // All ASCII characters allowed, except for {}
    BOOST_TEST(
        format_sql(opts, "{:abcdefghijklmnopqrstuvwxyz}", echo_spec{}) == "abcdefghijklmnopqrstuvwxyz"
    );
    BOOST_TEST(
        format_sql(opts, "{0:ABCDEFGHIJKLMNOPQRSTUVWXYZ}", echo_spec{}) == "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    );
    BOOST_TEST(format_sql(opts, "{name:0abAZ12 3456789}", named_args) == "0abAZ12 3456789");
    BOOST_TEST(
        format_sql(opts, "{:!\"#$%&'()*+,-./:;<=>?@[]\\^_`|~}", echo_spec{}) ==
        "!\"#$%&'()*+,-./:;<=>?@[]\\^_`|~"
    );
}

// Returning something != end in parse is an error
BOOST_AUTO_TEST_CASE(parse_error)
{
    BOOST_TEST(
        format_single_error("{:abc}", parse_consume<0>()) == client_errc::format_string_invalid_specifier
    );
    BOOST_TEST(
        format_single_error("{:abc}", parse_consume<1>()) == client_errc::format_string_invalid_specifier
    );
    BOOST_TEST(
        format_single_error("{:abc}", parse_consume<2>()) == client_errc::format_string_invalid_specifier
    );
    BOOST_TEST(format_single_error("{:abc}", parse_consume<3>()) == error_code());
}

// Formatters may accept values by const/non-const lvalue reference,
// and constness in format_sql is correctly propagated to them.
// This is used internally by sequence(), because some ranges are non-const.
// Not part of the public API - only const formatters are exposed
BOOST_AUTO_TEST_CASE(format_only_nonconst)
{
    // This would be the case of filter_view
    only_nonconst_lvalues lval{42};
    BOOST_TEST(format_sql(opts, "{}", lval) == "only_nonconst_lvalues(42)");
    BOOST_TEST(format_sql(opts, "{}", only_nonconst_lvalues{80}) == "only_nonconst_lvalues(80)");
}

BOOST_AUTO_TEST_CASE(format_const_nonconst)
{
    // Can appear in generic code
    const_and_nonconst lval{};
    const const_and_nonconst clval{};
    BOOST_TEST(format_sql(opts, "{}", lval) == "nonconst");
    BOOST_TEST(format_sql(opts, "{}", clval) == "const");
    BOOST_TEST(format_sql(opts, "{}", const_and_nonconst{}) == "nonconst");
    BOOST_TEST(format_sql(opts, "{}", std::move(clval)) == "const");
}

// format can call add_error to report problems
BOOST_AUTO_TEST_CASE(format_add_error)
{
    BOOST_TEST(format_single_error("SELECT {};", maybe_error{true}) == client_errc::unformattable_value);
    BOOST_TEST(format_sql(opts, "SELECT {};", maybe_error{false}) == "SELECT 0;");
}

// Spotcheck on a realistic type
BOOST_AUTO_TEST_CASE(spotcheck)
{
    BOOST_TEST(format_sql(opts, "{}", custom::condition{"myfield", 42}) == "`myfield`=42");
    BOOST_TEST(format_sql(opts, "{:s}", custom::condition{"myfield", 42}) == "`myfield` = 42");
    BOOST_TEST(
        format_single_error("{:i}", custom::condition{"myfield", 42}) ==
        client_errc::format_string_invalid_specifier
    );
}

BOOST_AUTO_TEST_SUITE_END()