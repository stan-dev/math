//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Antony Polukhin, 2011-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

#include "escape_struct.hpp"

#include <boost/core/lightweight_test.hpp>

#include <vector>

using namespace boost;

#if defined(BOOST_NO_STRINGSTREAM) || defined(BOOST_NO_STD_WSTRING)
#define BOOST_LCAST_NO_WCHAR_T
#endif

template <class T>
void do_test_on_empty_input(T& v)
{
    BOOST_TEST_THROWS(lexical_cast<int>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<float>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<double>(v), bad_lexical_cast);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    BOOST_TEST_THROWS(lexical_cast<long double>(v), bad_lexical_cast);
#endif
    BOOST_TEST_THROWS(lexical_cast<unsigned int>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<unsigned short>(v), bad_lexical_cast);
#if defined(BOOST_HAS_LONG_LONG)
    BOOST_TEST_THROWS(lexical_cast<boost::ulong_long_type>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<boost::long_long_type>(v), bad_lexical_cast);
#elif defined(BOOST_HAS_MS_INT64)
    BOOST_TEST_THROWS(lexical_cast<unsigned __int64>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<__int64>(v), bad_lexical_cast);
#endif
}

void test_empty_iterator_range()
{

    boost::iterator_range<char*> v;
    do_test_on_empty_input(v);
    BOOST_TEST_EQ(lexical_cast<std::string>(v), std::string());
    BOOST_TEST_THROWS(lexical_cast<char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<unsigned char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<signed char>(v), bad_lexical_cast);

    boost::iterator_range<const char*> cv;
    do_test_on_empty_input(cv);
    BOOST_TEST_EQ(lexical_cast<std::string>(cv), std::string());
    BOOST_TEST_THROWS(lexical_cast<char>(cv), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<unsigned char>(cv), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<signed char>(cv), bad_lexical_cast);

    const boost::iterator_range<const char*> ccv;
    do_test_on_empty_input(ccv);
    BOOST_TEST_EQ(lexical_cast<std::string>(ccv), std::string());
    BOOST_TEST_THROWS(lexical_cast<char>(ccv), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<unsigned char>(ccv), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<signed char>(ccv), bad_lexical_cast);
}

void test_empty_string()
{
    std::string v;
    do_test_on_empty_input(v);
    BOOST_TEST_THROWS(lexical_cast<char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<unsigned char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<signed char>(v), bad_lexical_cast);

#ifndef BOOST_LCAST_NO_WCHAR_T
    std::wstring vw;
    do_test_on_empty_input(vw);
    BOOST_TEST_THROWS(lexical_cast<wchar_t>(vw), bad_lexical_cast);
#endif

// Currently, no compiler and STL library fully support char16_t and char32_t
//#ifndef BOOST_NO_CXX11_CHAR16_T
//    std::basic_string<char16_t> v16w;
//    do_test_on_empty_input(v16w);
//    BOOST_TEST_THROWS(lexical_cast<char16_t>(v16w), bad_lexical_cast);
//#endif
//#ifndef BOOST_NO_CXX11_CHAR32_T
//    std::basic_string<char32_t> v32w;
//    do_test_on_empty_input(v32w);
//    BOOST_TEST_THROWS(lexical_cast<char32_t>(v32w), bad_lexical_cast);
//#endif
}

void test_empty_user_class()
{
    EscapeStruct v("");
    do_test_on_empty_input(v);
    BOOST_TEST_THROWS(lexical_cast<char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<unsigned char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<signed char>(v), bad_lexical_cast);
}

namespace std {
inline std::ostream & operator<<(std::ostream & out, const std::vector<long> & v)
{
    std::ostream_iterator<long> it(out);
    std::copy(v.begin(), v.end(), it);
    BOOST_TEST(out);
    return out;
}
}

void test_empty_vector()
{
    std::vector<long> v;
    do_test_on_empty_input(v);
    BOOST_TEST_THROWS(lexical_cast<char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<unsigned char>(v), bad_lexical_cast);
    BOOST_TEST_THROWS(lexical_cast<signed char>(v), bad_lexical_cast);
}


struct my_string {
    friend std::ostream &operator<<(std::ostream& sout, my_string const&/* st*/) {
            return sout << "";
    }
};

void test_empty_zero_terminated_string()
{
    my_string st;
    BOOST_TEST_EQ(boost::lexical_cast<std::string>(st), std::string());;
}


int main()
{
    test_empty_iterator_range();
    test_empty_string();
    test_empty_user_class();
    test_empty_vector();
    test_empty_zero_terminated_string();

    return boost::report_errors();
}
