//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Antony Polukhin, 2012-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).


#include <boost/lexical_cast.hpp>

#include <boost/utility/string_view.hpp>

#ifndef BOOST_NO_CXX17_HDR_STRING_VIEW
#include <string_view>
#endif

#include <boost/core/lightweight_test.hpp>


template <class StringView>
void test_string_view_conversion() {
    using boost::lexical_cast;

    StringView sw = "1";
    BOOST_TEST_EQ(lexical_cast<std::string>(sw), "1");
    BOOST_TEST_EQ(lexical_cast<int>(sw), 1);

    sw = StringView("a\0b", 3);
    BOOST_TEST_EQ(lexical_cast<std::string>(sw), std::string("a\0b", 3));

    sw = StringView("a\0b", 4);
    BOOST_TEST_EQ(lexical_cast<std::string>(sw), std::string("a\0b", 4));

    sw = StringView("\0a\0", 3);
    BOOST_TEST_EQ(lexical_cast<std::string>(sw), std::string("\0a\0", 3));
}

int main() {
    test_string_view_conversion<boost::string_view>();

#ifndef BOOST_NO_CXX17_HDR_STRING_VIEW
    test_string_view_conversion<std::string_view>();
#endif

    return boost::report_errors();
}
