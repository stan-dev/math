// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ extended_json_example
// This header includes a type called json::value that acts as a
// Javascript-like polymorphic value type.
#include "json.hpp"

#include <boost/parser/parser.hpp>

#include <fstream>
#include <vector>
#include <climits>


namespace json {

    namespace bp = ::boost::parser;
    using namespace bp::literals;

    // The JSON spec imposes a limit on how deeply JSON data structures are
    // allowed to nest.  This exception is thrown when that limit is exceeded
    // during the parse.
    template<typename Iter>
    struct excessive_nesting : std::runtime_error
    {
        excessive_nesting(Iter it) :
            runtime_error("excessive_nesting"),
            iter(it)
        {}
        Iter iter;
    };


    // The only globals we need to parse JSON are: "How many data structures
    // deep are we?", and "What is the limit of open data structures
    // allowed?".
    struct global_state
    {
        int recursive_open_count = 0;
        int max_recursive_open_count = 0;
    };

    // When matching paired UTF-16 surrogates, we need to track a bit of state
    // between matching the first and second UTF-16 code units: namely, the
    // value of the first code unit.
    struct double_escape_locals
    {
        int first_surrogate = 0;
    };


    // Here are all the rules declared.  I've given them names that are
    // end-user friendly, so that if there is a parse error, you get a message
    // like "expected four hexadecimal digits here:", instead of "expected
    // hex_4 here:".

    bp::rule<class ws> const ws = "whitespace";

    bp::rule<class string_char, uint32_t> const string_char =
        "code point (code points <= U+001F must be escaped)";
    bp::rule<class four_hex_digits, uint32_t> const hex_4 =
        "four hexadecimal digits";
    bp::rule<class escape_seq, uint32_t> const escape_seq =
        "\\uXXXX hexadecimal escape sequence";
    bp::rule<class escape_double_seq, uint32_t, double_escape_locals> const
        escape_double_seq = "\\uXXXX hexadecimal escape sequence";
    bp::rule<class single_escaped_char, uint32_t> const single_escaped_char =
        "'\"', '\\', '/', 'b', 'f', 'n', 'r', or 't'";

    bp::rule<class null, value> const null = "null";
    bp::rule<class string, std::string> const string = "string";
    bp::rule<class number, double> const number = "number";
    bp::rule<class object_element, boost::parser::tuple<std::string, value>> const
        object_element = "object-element";
    bp::rule<class object_tag, value> const object_p = "object";
    bp::rule<class array_tag, value> const array_p = "array";

    bp::rule<class value_tag, value> const value_p = "value";



    // JSON limits whitespace to just these four characters.
    auto const ws_def = '\x09'_l | '\x0a' | '\x0d' | '\x20';

    // Since our json object representation, json::value, is polymorphic, and
    // since its default-constructed state represents the JSON value "null",
    // we need to tell a json::value that it is an object (similar to a map)
    // before we start inserting values into it.  That's why we need
    // object_init.
    auto object_init = [](auto & ctx) {
        auto & globals = _globals(ctx);
        if (globals.max_recursive_open_count < ++globals.recursive_open_count)
            throw excessive_nesting(_where(ctx).begin());
        _val(ctx) = object();
    };

    // We need object_insert because we can't just insert into the json::value
    // itself.  The json::value does not have an insert() member, because if
    // it is currently holding a number, that makes no sense.  So, for a
    // json::value x, we need to call get<object>(x) to get the object
    // interface.
    auto object_insert = [](auto & ctx) {
        value & v = _val(ctx);
        get<object>(v).insert(std::make_pair(
            std::move(_attr(ctx))[0_c], std::move(_attr(ctx)[1_c])));
    };

    // These are the array analogues of the object semantic actions above.
    auto array_init = [](auto & ctx) {
        auto & globals = _globals(ctx);
        if (globals.max_recursive_open_count < ++globals.recursive_open_count)
            throw excessive_nesting(_where(ctx).begin());
        _val(ctx) = array();
    };
    auto array_append = [](auto & ctx) {
        value & v = _val(ctx);
        get<array>(v).push_back(std::move(_attr(ctx)));
    };

    // escape_double_seq is used to match pairs of UTF-16 surrogates that form
    // a single code point.  So, after matching one UTF-16 code unit c, we
    // only want to keep going if c is a lead/high surrogate.
    auto first_hex_escape = [](auto & ctx) {
        auto & locals = _locals(ctx);
        uint32_t const cu = _attr(ctx);
        if (!boost::parser::detail::text::high_surrogate(cu))
            _pass(ctx) = false; // Not a high surrogate; explicitly fail the parse.
        else
            locals.first_surrogate = cu; // Save this initial code unit for later.
    };
    // This is also used in escape_double_seq.  When we get to this action, we
    // know we've already matched a high surrogate, and so this one had better
    // be a low surrogate, or we have a (local) parse failure.
    auto second_hex_escape = [](auto & ctx) {
        auto & locals = _locals(ctx);
        uint32_t const cu = _attr(ctx);
        if (!boost::parser::detail::text::low_surrogate(cu)) {
            _pass(ctx) = false; // Not a low surrogate; explicitly fail the parse.
        } else {
            // Success!  Write to the rule's attribute the code point that the
            // first and second code points form.
            uint32_t const high_surrogate_min = 0xd800;
            uint32_t const low_surrogate_min = 0xdc00;
            uint32_t const surrogate_offset =
                0x10000 - (high_surrogate_min << 10) - low_surrogate_min;
            uint32_t const first_cu = locals.first_surrogate;
            _val(ctx) = (first_cu << 10) + cu + surrogate_offset;
        }
    };

    // This is the verbose form of declaration for the integer and unsigned
    // integer parsers int_parser and uint_parser.  In this case, we don't
    // want to use boost::parser::hex directly, since it has a variable number
    // of digits.  We want to match exactly 4 digits, and this is how we
    // declare a hexadecimal parser that matches exactly 4.
    bp::parser_interface<bp::uint_parser<uint32_t, 16, 4, 4>> const hex_4_def;

    // We use > here instead of >>, because once we see \u, we know that
    // exactly four hex digits must follow -- no other production rule starts
    // with \u.
    auto const escape_seq_def = "\\u" > hex_4;

    // This uses the actions above and the simpler rule escape_seq to find
    // matched UTF-16 surrogate pairs.
    auto const escape_double_seq_def =
        escape_seq[first_hex_escape] >> escape_seq[second_hex_escape];

    // This symbol table recognizes each character that can appear right after
    // an escaping backslash, and, if it finds one, produces the associated
    // code point as its attribute.
    bp::symbols<uint32_t> const single_escaped_char_def = {
        {"\"", 0x0022u},
        {"\\", 0x005cu},
        {"/", 0x002fu},
        {"b", 0x0008u},
        {"f", 0x000cu},
        {"n", 0x000au},
        {"r", 0x000du},
        {"t", 0x0009u}};

    // A string may be a matched UTF-16 escaped surrogate pair, a single
    // escaped UTF-16 code unit treated as a whole code point, a single
    // escaped character like \f, or any other code point outside the range
    // [0x0000u, 0x001fu].  Note that we had to put escape_double_seq before
    // escape_seq.  Otherwise, escape_seq would eat all the escape sequences
    // before escape_double_seq could try to match them.
    auto const string_char_def = escape_double_seq | escape_seq |
                                 ('\\'_l > single_escaped_char) |
                                 (bp::cp - bp::char_(0x0000u, 0x001fu));

    // If we see the special token null, treat that as a default-constructed
    // json::value.  Note that we could have done this with a semantic action,
    // but it is best to do everything you can without semantic actions;
    // they're a lot of code.
    auto const null_def = "null" >> bp::attr(value());

    auto const string_def = bp::lexeme['"' >> *(string_char - '"') > '"'];

    // Since the JSON format for numbers is not exactly what
    // boost::parser::double_ accepts (double_ accepts too much), we need to
    // parse a JSON number as a sequence of characters, and then pass the
    // result to double_ to actually get the numeric value.  This action does
    // that.  The parser uses boost::parser::raw to produce the subrange of
    // the input that covers the number as an attribute, which is used here.
    auto parse_double = [](auto & ctx) {
        auto const cp_range = _attr(ctx);
        auto cp_first = cp_range.begin();
        auto const cp_last = cp_range.end();

        auto const result = bp::prefix_parse(cp_first, cp_last, bp::double_);
        if (result) {
            _val(ctx) = *result;
        } else {
            // This would be more efficient if we used
            // boost::container::small_vector, or std::inplace_vector from
            // C++26.
            std::vector<char> chars(cp_first, cp_last);
            auto const chars_first = &*chars.begin();
            auto chars_last = chars_first + chars.size();
            _val(ctx) = std::strtod(chars_first, &chars_last);
        }
    };

    // As indicated above, we want to match the specific formats JSON allows,
    // and then re-parse the resulting matched range within the semantic
    // action.
    auto const number_def =
        bp::raw[bp::lexeme
                    [-bp::char_('-') >>
                     (bp::char_('1', '9') >> *bp::digit | bp::char_('0')) >>
                     -(bp::char_('.') >> +bp::digit) >>
                     -(bp::char_("eE") >> -bp::char_("+-") >> +bp::digit)]]
               [parse_double];

    // Note how, in the next three parsers, we turn off backtracking by using
    // > instead of >>, once we know that there is no backtracking alternative
    // that might match if we fail to match the next element.  This produces
    // much better error messages than if you always use >>.

    auto const object_element_def = string > ':' > value_p;

    auto const object_p_def = '{'_l[object_init] >>
                              -(object_element[object_insert] % ',') > '}';

    auto const array_p_def = '['_l[array_init] >>
                             -(value_p[array_append] % ',') > ']';

    // This is the top-level parser.
    auto const value_p_def =
        number | bp::bool_ | null | string | array_p | object_p;

    // Here, we define all the rules we've declared above, which also connects
    // each rule to its _def-suffixed parser.
    BOOST_PARSER_DEFINE_RULES(
        ws,
        hex_4,
        escape_seq,
        escape_double_seq,
        single_escaped_char,
        string_char,
        null,
        string,
        number,
        object_element,
        object_p,
        array_p,
        value_p);

    // json::parse() takes a string_view as input.  It takes an optional
    // callback to use for error reporting, which defaults to a no-op that
    // ignores all errors.  It also takes an optional max recursion depth
    // limit, which defaults to the one from the JSON spec, 512.
    std::optional<value> parse(
        std::string_view str,
        diagnostic_function errors_callback = diagnostic_function(),
        int max_recursion = 512)
    {
        // Turn the input range into a UTF-32 range, so that we can be sure
        // that we fall into the Unicode-aware parsing path inside parse()
        // below.
        auto const range = boost::parser::as_utf32(str);
        using iter_t = decltype(range.begin());

        if (max_recursion <= 0)
            max_recursion = INT_MAX;

        // Initialize our globals to the current depth (0), and the max depth
        // (max_recursion).
        global_state globals{0, max_recursion};
        bp::callback_error_handler error_handler(errors_callback);
        // Make a new parser that includes the globals and error handler.
        auto const parser = bp::with_error_handler(
            bp::with_globals(value_p, globals), error_handler);

        try {
            // Parse.  If no exception is thrown, due to: a failed expectation
            // (such as foo > bar, where foo matches the input, but then bar
            // cannot); or because the nesting depth is exceeded; we simply
            // return the result of the parse.  The result will contextually
            // convert to false if the parse failed.  Note that the
            // failed-expectation exception is caught internally, and used to
            // generate an error message.
            return bp::parse(range, parser, ws);
        } catch (excessive_nesting<iter_t> const & e) {
            // If we catch an excessive_nesting exception, just report it
            // and return an empty/failure result.
            if (errors_callback) {
                std::string const message = "error: Exceeded maximum number (" +
                                            std::to_string(max_recursion) +
                                            ") of open arrays and/or objects";
                std::stringstream ss;
                bp::write_formatted_message(
                    ss, "", range.begin(), e.iter, range.end(), message);
                errors_callback(ss.str());
            }
        }

        return {};
    }

}

std::string file_slurp(std::ifstream & ifs)
{
    std::string retval;
    while (ifs) {
        char const c = ifs.get();
        retval += c;
    }
    if (!retval.empty() && retval.back() == -1)
        retval.pop_back();
    return retval;
}

int main(int argc, char * argv[])
{
    if (argc < 2) {
        std::cerr << "A filename to parse is required.\n";
        exit(1);
    }

    std::ifstream ifs(argv[1]);
    if (!ifs) {
        std::cerr << "Unable to read file '" << argv[1] << "'.\n";
        exit(1);
    }

    // Read in the entire file.
    std::string const file_contents = file_slurp(ifs);
    // Parse the contents.  If there is an error, just stream it to cerr.
    auto json = json::parse(
        file_contents, [](std::string const & msg) { std::cerr << msg; });
    if (!json) {
        std::cerr << "Parse failure.\n";
        exit(1);
    }

    std::cout << "Parse successful; contents:\n" << *json << "\n";

    return 0;
}
//]
