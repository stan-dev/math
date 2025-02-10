// Copyright (C) 2018 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


using namespace boost::parser;

auto make_13 = [](auto & context) { return 13; };


auto const first_param_to_val = [](auto & context) {
    using namespace boost::parser::literals;
    _val(context) = (int)get(_params(context), 0_c);
};
auto const second_param_to_val = [](auto & context) {
    using namespace boost::parser::literals;
    _val(context) = get(_params(context), 1_c);
};
constexpr rule<struct action_param_tag, int> action_param = "abc or def";
constexpr auto action_param_def =
    string("abc")[first_param_to_val] | string("def")[second_param_to_val];
BOOST_PARSER_DEFINE_RULES(action_param);

int main()
{
    {
        std::string const str = "abc";
        auto const result = parse(str, action_param.with(15.0, make_13));
        BOOST_TEST(result);
        BOOST_TEST(*result == 15);
    }
    {
        std::string const str = "def";
        auto const result = parse(str, action_param.with(15.0, make_13));
        BOOST_TEST(result);
        BOOST_TEST(*result == 13);
    }

    return boost::report_errors();
}
