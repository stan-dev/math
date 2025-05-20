// Copyright (C) 2024 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//[ parsing_into_a_struct_example
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

int main()
{
    std::cout << "Enter employee record. ";
    std::string input;
    std::getline(std::cin, input);

    auto quoted_string = bp::lexeme['"' >> +(bp::char_ - '"') >> '"'];
    auto employee_p = bp::lit("employee")
        >> '{'
        >> bp::int_ >> ','
        >> quoted_string >> ','
        >> quoted_string >> ','
        >> bp::double_
        >> '}';

    employee record;
    auto const result = bp::parse(input, employee_p, bp::ws, record);

    if (result) {
        std::cout << "You entered:\nage:      " << record.age
                  << "\nsurname:  " << record.surname
                  << "\nforename: " << record.forename
                  << "\nsalary  : " << record.salary << "\n";
    } else {
        std::cout << "Parse failure.\n";
    }
}
//]
