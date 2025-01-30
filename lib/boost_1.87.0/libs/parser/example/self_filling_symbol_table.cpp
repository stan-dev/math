// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ self_filling_symbol_table_example
#include <boost/parser/parser.hpp>

#include <iostream>
#include <string>


namespace bp = boost::parser;

int main()
{
    //[ self_filling_symbol_table_table
    bp::symbols<int> const symbols = {{"c", 8}};
    assert(parse("c", symbols));
    //]

    //[ self_filling_symbol_table_action
    auto const add_symbol = [&symbols](auto & ctx) {
        using namespace bp::literals;
        // symbols::insert() requires a string, not a single character.
        char chars[2] = {_attr(ctx)[0_c], 0};
        symbols.insert(ctx, chars, _attr(ctx)[1_c]);
    };
    //]
    //[ self_filling_symbol_table_parser
    auto const parser = (bp::char_ >> bp::int_)[add_symbol] >> symbols;
    //]

    //[ self_filling_symbol_table_parse
    auto const result = parse("X 9 X", parser, bp::ws);
    assert(result && *result == 9);
    //]
    (void)result;

    //[ self_filling_symbol_table_after_parse
    assert(!parse("X", symbols));
    //]
}
//]
