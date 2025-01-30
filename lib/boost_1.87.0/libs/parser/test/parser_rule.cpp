// Copyright (C) 2018 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;

constexpr rule<struct flat_rule_tag> flat_rule = "flat_rule";
constexpr auto flat_rule_def = string("abc") | string("def");
BOOST_PARSER_DEFINE_RULES(flat_rule);

constexpr rule<struct recursive_rule_tag> recursive_rule = "recursive_rule";
constexpr auto recursive_rule_def = string("abc") >> -('a' >> recursive_rule);
BOOST_PARSER_DEFINE_RULES(recursive_rule);

constexpr rule<struct flat_string_rule_tag, std::string> flat_string_rule =
    "flat_string_rule";
constexpr auto flat_string_rule_def = string("abc") | string("def");
BOOST_PARSER_DEFINE_RULES(flat_string_rule);

constexpr callback_rule<struct recursive_string_rule_tag, std::string>
    recursive_string_rule = "recursive_string_rule";
auto append_string = [](auto & ctx) {
    auto & val = _val(ctx);
    auto & attr = _attr(ctx);
    val.insert(val.end(), attr.begin(), attr.end());
};
constexpr auto recursive_string_rule_def = string("abc")[append_string] >>
                                           -('a' >> recursive_string_rule);
BOOST_PARSER_DEFINE_RULES(recursive_string_rule);

constexpr rule<struct flat_vector_rule_tag, std::vector<char>>
    flat_vector_rule = "flat_vector_rule";
constexpr auto flat_vector_rule_def = string("abc") | string("def");
BOOST_PARSER_DEFINE_RULES(flat_vector_rule);

constexpr callback_rule<struct callback_vector_rule_tag, std::vector<char>>
    callback_vector_rule = "callback_vector_rule";
constexpr auto callback_vector_rule_def = string("abc") | string("def");
BOOST_PARSER_DEFINE_RULES(callback_vector_rule);

constexpr callback_rule<struct callback_void_rule_tag> callback_void_rule =
    "callback_void_rule";
constexpr auto callback_void_rule_def = string("abc") | string("def");
BOOST_PARSER_DEFINE_RULES(callback_void_rule);

struct callback_vector_rule_tag
{};
struct callback_void_rule_tag
{};

struct callbacks_t
{
    void operator()(callback_vector_rule_tag, std::vector<char> && vec) const
    {
        all_results.push_back(std::move(vec));
    }
    void operator()(callback_void_rule_tag) const
    {
        void_callback_called = true;
    }

    mutable std::vector<std::vector<char>> all_results;
    mutable bool void_callback_called = false;
};

struct recursive_strings_rule_tag
{};
constexpr callback_rule<recursive_strings_rule_tag, std::vector<std::string>>
    recursive_strings_rule = "recursive_strings_rule";
auto push_back = [](auto & ctx) { _val(ctx).push_back(std::move(_attr(ctx))); };
constexpr auto recursive_strings_rule_def = string("abc")[push_back] >>
                                            -('a' >> recursive_strings_rule);
BOOST_PARSER_DEFINE_RULES(recursive_strings_rule);

template<typename Tag, typename Attribute>
struct recursive_rule_callbacks_t
{
    void operator()(Tag, Attribute && vec) const
    {
        all_results.push_back(std::move(vec));
    }

    mutable std::vector<Attribute> all_results;
};

struct recursive_string_rule_tag
{};

namespace more_about_rules_1 {
    namespace bp = boost::parser;

    bp::rule<struct value_tag> value =
        "an integer, or a list of integers in braces";

    auto const ints = '{' > (value % ',') > '}';
    auto const value_def = bp::int_ | ints;

    BOOST_PARSER_DEFINE_RULES(value);
}

namespace more_about_rules_2 {
    namespace bp = boost::parser;

    bp::rule<struct value_tag> value =
        "an integer, or a list of integers in braces";
    bp::rule<struct comma_values_tag> comma_values =
        "a comma-delimited list of integers";

    auto const ints = '{' > comma_values > '}';
    auto const value_def = bp::int_ | ints;
    auto const comma_values_def = (value % ',');

    BOOST_PARSER_DEFINE_RULES(value, comma_values);
}

namespace more_about_rules_3 {
    namespace bp = boost::parser;

    bp::rule<struct parens_tag> parens = "matched parentheses";

    auto const parens_def = (('(' >> parens) > ')') | bp::eps;

    BOOST_PARSER_DEFINE_RULES(parens);
}

namespace more_about_rules_4 {
    namespace bp = boost::parser;

    bp::rule<struct parens_tag> parens = "matched parentheses";

    auto const parens_def = (('(' >> parens) > ')') | bp::eps;

    BOOST_PARSER_DEFINE_RULES(parens);

    bool balanced_parens(std::string_view str);
}

namespace more_about_rules_4 {
    bool balanced_parens(std::string_view str)
    {
        namespace bp = boost::parser;
        return bp::parse(str, bp::omit[parens], bp::ws);
    }
}

// clang-format off
namespace param_example {
    //[ extended_param_yaml_example_rules
    namespace bp = boost::parser;

    // A type to represent the YAML parse context.
    enum class context {
        block_in,
        block_out,
        block_key,
        flow_in,
        flow_out,
        flow_key
    };

    // A YAML value; no need to fill it in for this example.
    struct value
    {
        // ...
    };

    // YAML [66], just stubbed in here.
    auto const s_separate_in_line = bp::eps;

    // YAML [137].
    bp::rule<struct c_flow_seq_tag, value> c_flow_sequence = "c-flow-sequence";
    // YAML [80].
    bp::rule<struct s_separate_tag> s_separate = "s-separate";
    // YAML [136].
    bp::rule<struct in_flow_tag, value> in_flow = "in-flow";
    // YAML [138]; just eps below.
    bp::rule<struct ns_s_flow_seq_entries_tag, value> ns_s_flow_seq_entries =
        "ns-s-flow-seq-entries";
    // YAML [81]; just eps below.
    bp::rule<struct s_separate_lines_tag> s_separate_lines = "s-separate-lines";

    // Parser for YAML [137].
    auto const c_flow_sequence_def =
        '[' >>
        -s_separate.with(bp::_p<0>, bp::_p<1>) >>
        -in_flow.with(bp::_p<0>, bp::_p<1>) >>
        ']';
    // Parser for YAML [80].
    auto const s_separate_def = bp::switch_(bp::_p<1>)
        (context::block_out, s_separate_lines.with(bp::_p<0>))
        (context::block_in, s_separate_lines.with(bp::_p<0>))
        (context::flow_out, s_separate_lines.with(bp::_p<0>))
        (context::flow_in, s_separate_lines.with(bp::_p<0>))
        (context::block_key, s_separate_in_line)
        (context::flow_key, s_separate_in_line);
    // Parser for YAML [136].
    auto const in_flow_def = bp::switch_(bp::_p<1>)
        (context::flow_out, ns_s_flow_seq_entries.with(bp::_p<0>, context::flow_in))
        (context::flow_in, ns_s_flow_seq_entries.with(bp::_p<0>, context::flow_in))
        (context::block_out, ns_s_flow_seq_entries.with(bp::_p<0>, context::flow_key))
        (context::flow_key, ns_s_flow_seq_entries.with(bp::_p<0>, context::flow_key));

    auto const ns_s_flow_seq_entries_def = bp::eps;
    auto const s_separate_lines_def = bp::eps;

    BOOST_PARSER_DEFINE_RULES(
        c_flow_sequence,
        s_separate,
        in_flow,
        ns_s_flow_seq_entries,
        s_separate_lines);
//]
}
// clang-format on

int main()
{

// no_attribute_rules
{
    {
        std::string const str = "xyz";
        BOOST_TEST(!parse(str, flat_rule));
        BOOST_TEST(!parse(str, recursive_rule));
    }
    {
        std::string const str = "def";
        bool const flat_result{parse(str, flat_rule)};
        BOOST_TEST(flat_result);
        BOOST_TEST(!parse(str, recursive_rule));
    }
    {
        std::string const str = "abc";
        BOOST_TEST(parse(str, flat_rule));
        BOOST_TEST(parse(str, recursive_rule));
    }
    {
        std::string const str = "abcaabc";
        BOOST_TEST(!parse(str, flat_rule));
        BOOST_TEST(parse(str, recursive_rule));
    }
    {
        std::string const str = "abcaabc";
        auto first = str.c_str();
        BOOST_TEST(prefix_parse(first, boost::parser::detail::text::null_sentinel, flat_rule));
        first = str.c_str();
        BOOST_TEST(prefix_parse(first, boost::parser::detail::text::null_sentinel, recursive_rule));
    }
}

// string_attribute_rules
{
    {
        std::string const str = "xyz";
        BOOST_TEST(!parse(str, flat_string_rule));
        BOOST_TEST(!parse(str, recursive_string_rule));
    }
    {
        std::string const str = "def";
        auto const flat_result = parse(str, flat_string_rule);
        BOOST_TEST(flat_result);
        BOOST_TEST(*flat_result == "def");
        BOOST_TEST(!parse(str, recursive_string_rule));
    }
    {
        std::string const str = "abc";
        BOOST_TEST(*parse(str, flat_string_rule) == "abc");
        BOOST_TEST(*parse(str, recursive_string_rule) == "abc");
    }
    {
        std::string const str = "abcaabc";
        BOOST_TEST(!parse(str, flat_string_rule));
        BOOST_TEST(parse(str, recursive_string_rule));
    }
    {
        std::string const str = "abcaabc";
        auto first = str.c_str();
        BOOST_TEST(
            *prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                flat_string_rule) == "abc");
        first = str.c_str();
        BOOST_TEST(
            *prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                recursive_string_rule) == "abcabc");
    }
}

// vector_attribute_rules
{
    {
        std::string const str = "xyz";
        std::vector<char> chars;
        BOOST_TEST(!parse(str, flat_vector_rule, chars));
    }
    {
        std::string const str = "def";
        std::vector<char> chars;
        BOOST_TEST(parse(str, flat_vector_rule, chars));
        BOOST_TEST(chars == std::vector<char>({'d', 'e', 'f'}));
    }
    {
        std::string const str = "abc";
        BOOST_TEST(
            *parse(str, flat_vector_rule) ==
            std::vector<char>({'a', 'b', 'c'}));
    }
    {
        std::string const str = "abcaabc";
        BOOST_TEST(!parse(str, flat_vector_rule));
    }
    {
        std::string const str = "abcaabc";
        auto first = str.c_str();
        BOOST_TEST(
            *prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                flat_vector_rule) == std::vector<char>({'a', 'b', 'c'}));
    }
    {
        std::string const str = "abcaabc";
        auto first = str.c_str();
        BOOST_TEST(
            callback_prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                flat_vector_rule,
                int{}) == true);
    }
}

// callback_rules
{
    {
        std::string const str = "xyz";
        callbacks_t callbacks;
        BOOST_TEST(!callback_parse(str, callback_vector_rule, callbacks));
        BOOST_TEST(callbacks.all_results.size() == 0u);
    }
    {
        std::string const str = "abc";
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(str, callback_vector_rule, callbacks));
        BOOST_TEST(callbacks.all_results.size() == 1u);
        BOOST_TEST(callbacks.all_results[0] == std::vector<char>({'a', 'b', 'c'}));
    }
    {
        std::string const str = "def";
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(str, callback_vector_rule, callbacks));
        BOOST_TEST(callbacks.all_results.size() == 1u);
        BOOST_TEST(callbacks.all_results[0] == std::vector<char>({'d', 'e', 'f'}));
    }

    {
        std::string const str = "xyz";
        callbacks_t callbacks;
        BOOST_TEST(!callback_parse(str, callback_void_rule, callbacks));
        BOOST_TEST(!callbacks.void_callback_called);
    }
    {
        std::string const str = "abc";
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(str, callback_void_rule, callbacks));
        BOOST_TEST(callbacks.void_callback_called);
    }
    {
        std::string const str = "def";
        callbacks_t callbacks;
        BOOST_TEST(callback_parse(str, callback_void_rule, callbacks));
        BOOST_TEST(callbacks.void_callback_called);
    }

    {
        std::string const str = "xyz";
        std::vector<std::vector<char>> all_results;
        auto callbacks =
            [&all_results](callback_vector_rule_tag, std::vector<char> && vec) {
                all_results.push_back(std::move(vec));
            };
        BOOST_TEST(!callback_parse(str, callback_vector_rule, callbacks));
        BOOST_TEST(all_results.size() == 0u);
    }
    {
        std::string const str = "abc";
        std::vector<std::vector<char>> all_results;
        auto callbacks =
            [&all_results](callback_vector_rule_tag, std::vector<char> && vec) {
                all_results.push_back(std::move(vec));
            };
        BOOST_TEST(callback_parse(str, callback_vector_rule, callbacks));
        BOOST_TEST(all_results.size() == 1u);
        BOOST_TEST(all_results[0] == std::vector<char>({'a', 'b', 'c'}));
    }
    {
        std::string const str = "def";
        std::vector<std::vector<char>> all_results;
        auto callbacks =
            [&all_results](callback_vector_rule_tag, std::vector<char> && vec) {
                all_results.push_back(std::move(vec));
            };
        BOOST_TEST(callback_parse(str, callback_vector_rule, callbacks));
        BOOST_TEST(all_results.size() == 1u);
        BOOST_TEST(all_results[0] == std::vector<char>({'d', 'e', 'f'}));
    }

    {
        std::string const str = "xyz";
        bool void_callback_called = false;
        auto callbacks = [&void_callback_called](callback_void_rule_tag) {
            void_callback_called = true;
        };
        BOOST_TEST(!callback_parse(str, callback_void_rule, callbacks));
        BOOST_TEST(!void_callback_called);
    }
    {
        std::string const str = "abc";
        bool void_callback_called = false;
        auto callbacks = [&void_callback_called](callback_void_rule_tag) {
            void_callback_called = true;
        };
        BOOST_TEST(callback_parse(str, callback_void_rule, callbacks));
        BOOST_TEST(void_callback_called);
    }
    {
        std::string const str = "def";
        bool void_callback_called = false;
        auto callbacks = [&void_callback_called](callback_void_rule_tag) {
            void_callback_called = true;
        };
        BOOST_TEST(callback_parse(str, callback_void_rule, callbacks));
        BOOST_TEST(void_callback_called);
    }
}

// callback_rules_normal_parse
{
    {
        std::string const str = "xyz";
        std::vector<char> chars;
        BOOST_TEST(!parse(str, callback_vector_rule, chars));
    }
    {
        std::string const str = "abc";
        std::vector<char> chars;
        BOOST_TEST(parse(str, callback_vector_rule, chars));
        BOOST_TEST(chars == std::vector<char>({'a', 'b', 'c'}));
    }
    {
        std::string const str = "def";
        std::vector<char> chars;
        BOOST_TEST(parse(str, callback_vector_rule, chars));
        BOOST_TEST(chars == std::vector<char>({'d', 'e', 'f'}));
    }

    {
        std::string const str = "def";
        BOOST_TEST(parse(str, callback_void_rule));
    }

    {
        std::string const str = "xyz";
        auto const chars = parse(str, callback_vector_rule);
        BOOST_TEST(!chars);
    }
    {
        std::string const str = "abc";
        auto const chars = parse(str, callback_vector_rule);
        BOOST_TEST(chars);
        BOOST_TEST(chars == std::vector<char>({'a', 'b', 'c'}));
    }
    {
        std::string const str = "def";
        auto const chars = parse(str, callback_vector_rule);
        BOOST_TEST(chars);
        BOOST_TEST(chars == std::vector<char>({'d', 'e', 'f'}));
    }

    {
        std::string const str = "xyz";
        BOOST_TEST(!parse(str, callback_vector_rule));
    }
    {
        std::string const str = "abc";
        BOOST_TEST(parse(str, callback_vector_rule));
    }
    {
        std::string const str = "def";
        BOOST_TEST(parse(str, callback_vector_rule));
    }
}

////////////////////////////////////////////////////////////////////////////////
// More recursive rules

// container_populating_recursive_rule
{
    {
        std::string const str = "xyz";
        BOOST_TEST(!parse(str, recursive_strings_rule));
    }
    {
        std::string const str = "abc";
        BOOST_TEST(
            *parse(str, recursive_strings_rule) ==
            std::vector<std::string>({"abc"}));
    }
    {
        std::string const str = "abcaabc";
        BOOST_TEST(parse(str, recursive_strings_rule));
    }
    {
        std::string const str = "abcaabc";
        auto first = str.c_str();
        BOOST_TEST(
            *prefix_parse(
                first,
                boost::parser::detail::text::null_sentinel,
                recursive_strings_rule) ==
            std::vector<std::string>({"abc", "abc"}));
    }
}

// container_populating_recursive_cb_rule
{
    using callbacks_type = recursive_rule_callbacks_t<
        recursive_strings_rule_tag,
        std::vector<std::string>>;

    {
        std::string const str = "xyz";
        callbacks_type callbacks;
        BOOST_TEST(!callback_parse(str, recursive_strings_rule, callbacks));
    }
    {
        std::string const str = "abc";
        callbacks_type callbacks;
        BOOST_TEST(callback_parse(str, recursive_strings_rule, callbacks));
        BOOST_TEST(callbacks.all_results.size() == 1u);
        BOOST_TEST(callbacks.all_results[0] == std::vector<std::string>({"abc"}));
    }
    {
        std::string const str = "abcaabc";
        callbacks_type callbacks;
        BOOST_TEST(callback_parse(str, recursive_strings_rule, callbacks));
    }
    {
        std::string const str = "abcaabc";
        auto first = str.c_str();
        callbacks_type callbacks;
        BOOST_TEST(callback_prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            recursive_strings_rule,
            callbacks));
        BOOST_TEST(callbacks.all_results.size() == 1u);
        BOOST_TEST(
            callbacks.all_results[0] ==
            std::vector<std::string>({"abc", "abc"}));
    }
}

// string_recursive_cb_rule
{
    using callbacks_type =
        recursive_rule_callbacks_t<recursive_string_rule_tag, std::string>;

    {
        std::string const str = "xyz";
        callbacks_type callbacks;
        BOOST_TEST(!callback_parse(str, recursive_string_rule, callbacks));
    }
    {
        std::string const str = "def";
        callbacks_type callbacks;
        BOOST_TEST(!callback_parse(str, recursive_string_rule, callbacks));
    }
    {
        std::string const str = "abc";
        callbacks_type callbacks;
        BOOST_TEST(callback_parse(str, recursive_string_rule, callbacks));
        BOOST_TEST(callbacks.all_results.size() == 1u);
        BOOST_TEST(callbacks.all_results[0] == "abc");
    }
    {
        std::string const str = "abcaabc";
        callbacks_type callbacks;
        BOOST_TEST(callback_parse(str, recursive_string_rule, callbacks));
    }
    {
        std::string const str = "abcaabc";
        auto first = str.c_str();
        callbacks_type callbacks;
        BOOST_TEST(callback_prefix_parse(
            first,
            boost::parser::detail::text::null_sentinel,
            recursive_string_rule,
            callbacks));
        BOOST_TEST(callbacks.all_results.size() == 1u);
        BOOST_TEST(callbacks.all_results[0] == "abcabc");
    }
}

// doc_exaparensles
{
    {
        using namespace more_about_rules_1;

        bp::parse("{ 4, 5 a", value, bp::ws);
        bp::parse("{ }", value, bp::ws);
    }
    {
        using namespace more_about_rules_2;

        bp::parse("{ }", value, bp::ws);
    }
    {
        using namespace more_about_rules_3;

        BOOST_TEST(bp::parse("(((())))", parens, bp::ws));
        BOOST_TEST(!bp::parse("(((()))", parens, bp::ws));
    }
    {
        using namespace more_about_rules_4;
        assert(balanced_parens("(())"));
    }
}

// extended_param_example
{
    using namespace param_example;

    //[ extended_param_yaml_example_use
    auto const test_parser = c_flow_sequence.with(4, context::block_out);
    auto result = bp::parse("[]", test_parser);
    assert(result);
    //]
    (void)result;
}

return boost::report_errors();
}
