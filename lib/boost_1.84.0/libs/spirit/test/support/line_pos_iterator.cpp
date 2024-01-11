//  Copyright (c) 2014 Tomoki Imai
//  Copyright (c) 2023 Nikita Knizev
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying 
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/spirit/include/support_line_pos_iterator.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/assign.hpp>
#include <boost/static_assert.hpp>
#include <iterator>
#include <string>
#include <vector>

struct validation {
    validation(size_t line, size_t column)
    : line(line), column(column), inside_newline(true) {
    }

    validation(size_t line, size_t column, std::string current)
    : line(line), column(column), current(current), inside_newline(false) {
    }

    size_t line;
    size_t column;
    std::string current;
    bool inside_newline;
};

typedef std::vector<validation> validations;

struct string_forward_iterator : public boost::iterator_adaptor<
    string_forward_iterator      // Derived
  , std::string::const_iterator  // Base
  , boost::use_default           // Value
  , boost::forward_traversal_tag // CategoryOrTraversal
  >
{
    explicit string_forward_iterator(std::string::const_iterator it)
      : string_forward_iterator::iterator_adaptor_(it) {}
};
BOOST_STATIC_ASSERT(boost::is_same<std::iterator_traits<string_forward_iterator>::iterator_category, std::forward_iterator_tag>::value);

template <typename Iterator>
void test(Iterator first, Iterator last, validations const& validations, bool singlechar_newline) {
    typedef boost::spirit::line_pos_iterator<Iterator> pos_iterator_t;

    pos_iterator_t const input_begin(first);
    pos_iterator_t const input_end(last);
    pos_iterator_t position(input_begin);
    validations::const_iterator expected = validations.begin();

    for (; position != input_end && expected != validations.end(); ++position, ++expected) {
        if (expected->inside_newline && singlechar_newline)
            if (++expected == validations.end())
                break;
        boost::iterator_range<pos_iterator_t> const range = get_current_line(input_begin, position, input_end);
        std::string const current(range.begin(), range.end());

        BOOST_TEST_EQ(expected->line, get_line(position));
        BOOST_TEST_EQ(expected->column, get_column(input_begin, position));
        BOOST_TEST_EQ(expected->current, current);
    }
    if (expected != validations.end() && expected->inside_newline && singlechar_newline)
        ++expected;

    BOOST_TEST(position == input_end);
    BOOST_TEST(expected == validations.end());
}

void test(std::string const& input, validations const& validations, bool singlechar_newline) {
    test(input.begin(), input.end(), validations, singlechar_newline);
    test(string_forward_iterator(input.begin()), string_forward_iterator(input.end()),
         validations, singlechar_newline);
}

void test_only(std::string const& line_break) {
    std::string const input = line_break + line_break + line_break;
    validations const validations = boost::assign::list_of<validation>
        (1,1,"")(2,1)
        (2,1,"")(3,1)
        (3,1,"")(4,1);
    test(input, validations, line_break.size() == 1);
}

void test_in_text(std::string const& line_break) {
    std::string const input = "foo" + line_break + "bar" + line_break + "git";
    validations const validations = boost::assign::list_of<validation>
        (1,1,"foo")(1,2,"foo")(1,3,"foo")(1,4,"foo")(2,1)
        (2,1,"bar")(2,2,"bar")(2,3,"bar")(2,4,"bar")(3,1)
        (3,1,"git")(3,2,"git")(3,3,"git");
    test(input, validations, line_break.size() == 1);
}

void test_at_begin(std::string const& line_break) {
    std::string const input = line_break + "bar" + line_break + "git";
    validations const validations = boost::assign::list_of<validation>
        (1,1,"")(2,1)
        (2,1,"bar")(2,2,"bar")(2,3,"bar")(2,4,"bar")(3,1)
        (3,1,"git")(3,2,"git")(3,3,"git");
    test(input, validations, line_break.size() == 1);
}

void test_at_end(std::string const& line_break) {
    std::string const input = "foo" + line_break + "bar" + line_break;
    validations const validations = boost::assign::list_of<validation>
        (1,1,"foo")(1,2,"foo")(1,3,"foo")(1,4,"foo")(2,1)
        (2,1,"bar")(2,2,"bar")(2,3,"bar")(2,4,"bar")(3,1);
    test(input, validations, line_break.size() == 1);
}

void test_mixed() {
    std::string const input = "\n\n\r\r\n\r";
    validations const validations = boost::assign::list_of<validation>
        (1,1,"")      // \n
        (2,1,"")(3,1) // \n\r
        (3,1,"")(4,1) // \r\n
        (4,1,"");     // \r
    test(input, validations, false);
}

void test(std::string const& line_break)
{
    test_only(line_break);
    test_in_text(line_break);
    test_at_begin(line_break);
    test_at_end(line_break);
}

int main()
{
    test("\n");
    test("\r\n");
    test("\r");
    test("\n\r");
    test_mixed();

    return boost::report_errors();
}
