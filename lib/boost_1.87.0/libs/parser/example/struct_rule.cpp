// Copyright (C) 2024 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ struct_rule_example
#include <boost/parser/parser.hpp>

#include <iostream>
#include <string>


struct employee
{
    int age;
    std::string surname;
    std::string forename;
    double salary;
};

namespace bp = boost::parser;

bp::rule<struct quoted_string, std::string> quoted_string = "quoted name";
bp::rule<struct employee_p, employee> employee_p = "employee";

auto quoted_string_def = bp::lexeme['"' >> +(bp::char_ - '"') >> '"'];
auto employee_p_def = bp::lit("employee")
    >> '{'
    >> bp::int_ >> ','
    >> quoted_string >> ','
    >> quoted_string >> ','
    >> bp::double_
    >> '}';

BOOST_PARSER_DEFINE_RULES(quoted_string, employee_p);

int main()
{
    std::cout << "Enter employee record. ";
    std::string input;
    std::getline(std::cin, input);

    static_assert(std::is_aggregate_v<std::decay_t<employee &>>);

    auto const result = bp::parse(input, employee_p, bp::ws);

    if (result) {
        std::cout << "You entered:\nage:      " << result->age
                  << "\nsurname:  " << result->surname
                  << "\nforename: " << result->forename
                  << "\nsalary  : " << result->salary << "\n";
    } else {
        std::cout << "Parse failure.\n";
    }
}
//]
