// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ extended_callback_parsing_json_example
#include <boost/parser/parser.hpp>
#include <boost/parser/transcode_view.hpp>

#include <fstream>
#include <vector>
#include <climits>


namespace json {

    namespace bp = ::boost::parser;
    using namespace bp::literals;

    template<typename Iter>
    struct excessive_nesting : std::runtime_error
    {
        excessive_nesting(Iter it) :
            runtime_error("excessive_nesting"), iter(it)
        {}
        Iter iter;
    };


    struct global_state
    {
        int recursive_open_count = 0;
        int max_recursive_open_count = 0;
    };

    struct double_escape_locals
    {
        int first_surrogate = 0;
    };


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

    bp::callback_rule<class null_tag> const null = "null";

    // Since we don't create polymorphic values in this parse, we need to be
    // able to report that we parsed a bool, so we need a callback rule for
    // this.
    bp::callback_rule<class bool_tag, bool> const bool_p = "boolean";

    bp::callback_rule<class string_tag, std::string> const string = "string";
    bp::callback_rule<class number_tag, double> const number = "number";

    // object_element is broken up into the key (object_element_key) and the
    // whole thing (object_element).  This was done because the value after
    // the ':' may have many parts.  It may be an array, for example.  This
    // implies that we need to report that we have the string part of the
    // object-element, and that the rest -- the value -- is coming.
    bp::callback_rule<class object_element_key_tag, std::string> const
        object_element_key = "string";
    bp::rule<class object_element_tag> const object_element = "object-element";

    // object gets broken up too, to enable the reporting of the beginning and
    // end of the object when '{' or '}' is parsed, respectively.  The same
    // thing is done for array, below.
    bp::callback_rule<class object_open_tag> const object_open = "'{'";
    bp::callback_rule<class object_close_tag> const object_close = "'}'";
    bp::rule<class object_tag> const object = "object";

    bp::callback_rule<class array_open_tag> const array_open = "'['";
    bp::callback_rule<class array_close_tag> const array_close = "']'";
    bp::rule<class array_tag> const array = "array";

    // value no longer produces an attribute, and it has no callback either.
    // Each individual possible kind of value (string, array, etc.) gets
    // reported separately.
    bp::rule<class value_tag> const value = "value";


    // Since we use these tag types as function parameters in the callbacks,
    // they need to be complete types.
    class null_tag {};
    class bool_tag {};
    class string_tag {};
    class number_tag {};
    class object_element_key_tag {};
    class object_open_tag {};
    class object_close_tag {};
    class array_open_tag {};
    class array_close_tag {};


    auto const ws_def = '\x09'_l | '\x0a' | '\x0d' | '\x20';

    auto first_hex_escape = [](auto & ctx) {
        auto & locals = _locals(ctx);
        uint32_t const cu = _attr(ctx);
        if (!boost::parser::detail::text::high_surrogate(cu))
            _pass(ctx) = false;
        else
            locals.first_surrogate = cu;
    };
    auto second_hex_escape = [](auto & ctx) {
        auto & locals = _locals(ctx);
        uint32_t const cu = _attr(ctx);
        if (!boost::parser::detail::text::low_surrogate(cu)) {
            _pass(ctx) = false;
        } else {
            uint32_t const high_surrogate_min = 0xd800;
            uint32_t const low_surrogate_min = 0xdc00;
            uint32_t const surrogate_offset =
                0x10000 - (high_surrogate_min << 10) - low_surrogate_min;
            uint32_t const first_cu = locals.first_surrogate;
            _val(ctx) = (first_cu << 10) + cu + surrogate_offset;
        }
    };

    bp::parser_interface<bp::uint_parser<uint32_t, 16, 4, 4>> const hex_4_def;

    auto const escape_seq_def = "\\u" > hex_4;

    auto const escape_double_seq_def =
        escape_seq[first_hex_escape] >> escape_seq[second_hex_escape];

    bp::symbols<uint32_t> const single_escaped_char_def = {
        {"\"", 0x0022u},
        {"\\", 0x005cu},
        {"/", 0x002fu},
        {"b", 0x0008u},
        {"f", 0x000cu},
        {"n", 0x000au},
        {"r", 0x000du},
        {"t", 0x0009u}};

    auto const string_char_def = escape_double_seq | escape_seq |
                                 ('\\'_l > single_escaped_char) |
                                 (bp::cp - bp::char_(0x0000u, 0x001fu));

    auto const null_def = "null"_l;

    auto const bool_p_def = bp::bool_;

    auto const string_def = bp::lexeme['"' >> *(string_char - '"') > '"'];

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

    auto const number_def =
        bp::raw[bp::lexeme
                    [-bp::char_('-') >>
                     (bp::char_('1', '9') >> *bp::digit | bp::char_('0')) >>
                     -(bp::char_('.') >> +bp::digit) >>
                     -(bp::char_("eE") >> -bp::char_("+-") >> +bp::digit)]]
               [parse_double];

    // The object_element_key parser is exactly the same as the string parser.
    // Note that we did *not* use string here, though; we used string_def.  If
    // we had used string, its callback would have been called first, and
    // worse still, since it moves its attribute, the callback for
    // object_element_key would always report the empty string, because the
    // string callback would have consumed it first.
    auto const object_element_key_def = string_def;

    auto const object_element_def = object_element_key > ':' > value;

    // This is a very straightforward way to write object_def when we know we
    // don't care about attribute-generating (non-callback) parsing.  If we
    // wanted to support both modes in one parser definition, we could have
    // written:
    //    auto const object_open_def = eps;
    //    auto const object_close_def = eps;
    //    auto const object_def = '{' >> object_open >>
    //                             -(object_element % ',') >
    //                            '}' >> object_close;
    auto const object_open_def = '{'_l;
    auto const object_close_def = '}'_l;
    auto const object_def = object_open >>
                            -(object_element % ',') > object_close;

    auto const array_open_def = '['_l;
    auto const array_close_def = ']'_l;
    auto const array_def = array_open >> -(value % ',') > array_close;

    auto const value_def = number | bool_p | null | string | array | object;

    BOOST_PARSER_DEFINE_RULES(
        ws,
        hex_4,
        escape_seq,
        escape_double_seq,
        single_escaped_char,
        string_char,
        null,
        bool_p,
        string,
        number,
        object_element_key,
        object_element,
        object_open,
        object_close,
        object,
        array_open,
        array_close,
        array,
        value);

    // The parse function loses its attribute from the return type; now the
    // return type is just bool.
    template<typename Callbacks>
    bool parse(
        std::string_view str,
        std::string_view filename,
        Callbacks const & callbacks,
        int max_recursion = 512)
    {
        auto const range = boost::parser::as_utf32(str);
        using iter_t = decltype(range.begin());

        if (max_recursion <= 0)
            max_recursion = INT_MAX;

        global_state globals{0, max_recursion};
        // This is a different error handler from the json.cpp example, just
        // to show different options.
        bp::stream_error_handler error_handler(filename);
        auto const parser = bp::with_error_handler(
            bp::with_globals(value, globals), error_handler);

        try {
            // This is identical to the parse() call in json.cpp, except that
            // it is callback_parse() instead, and it takes the callbacks
            // parameter.
            return bp::callback_parse(range, parser, ws, callbacks);
        } catch (excessive_nesting<iter_t> const & e) {
            std::string const message = "error: Exceeded maximum number (" +
                                        std::to_string(max_recursion) +
                                        ") of open arrays and/or objects";
            bp::write_formatted_message(
                std::cout,
                filename,
                range.begin(),
                e.iter,
                range.end(),
                message);
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

// This is our callbacks-struct.  It has a callback for each of the kinds of
// callback rules in our parser.  If one were missing, you'd get a pretty
// nasty template instantiation error.  Note that these are all const members;
// callback_parse() takes the callbacks object by constant reference.
struct json_callbacks
{
    void operator()(json::null_tag) const { std::cout << "JSON null value\n"; }
    void operator()(json::bool_tag, bool b) const
    {
        indent();
        std::cout << "JSON bool " << (b ? "true" : "false") << "\n";
    }
    void operator()(json::string_tag, std::string s) const
    {
        indent();
        std::cout << "JSON string \"" << s << "\"\n";
    }
    void operator()(json::number_tag, double d) const
    {
        indent();
        std::cout << "JSON number " << d << "\n";
    }
    void operator()(json::object_element_key_tag, std::string key) const
    {
        indent();
        std::cout << "JSON object element with key \"" << key
                  << "\" and value...\n";
    }
    void operator()(json::object_open_tag) const
    {
        indent(1);
        std::cout << "Beginning of JSON object.\n";
    }
    void operator()(json::object_close_tag) const
    {
        indent(-1);
        std::cout << "End of JSON object.\n";
    }
    void operator()(json::array_open_tag) const
    {
        indent(1);
        std::cout << "Beginning of JSON array.\n";
    }
    void operator()(json::array_close_tag) const
    {
        indent(-1);
        std::cout << "End of JSON array.\n";
    }

    void indent(int level_bump = 0) const
    {
        if (level_bump < 0)
            indent_.resize(indent_.size() - 2);
        std::cout << indent_;
        if (0 < level_bump)
            indent_ += "  ";
    }
    mutable std::string indent_;
};

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

    std::string const file_contents = file_slurp(ifs);
    bool success = json::parse(file_contents, argv[1], json_callbacks{});
    if (success) {
        std::cout << "Parse successful!\n";
    } else {
        std::cerr << "Parse failure.\n";
        exit(1);
    }

    return 0;
}
//]
