//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Alexander Nasonov, 2007.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).
//
// This tests now must pass on vc8, because lexical_cast
// implementation has changed and it does not use stringstream for casts
// to integral types

#include <boost/lexical_cast.hpp>
#include <boost/cstdint.hpp>

#include <boost/core/lightweight_test.hpp>

#include <string>

#ifdef BOOST_MSVC
# pragma warning(disable: 4127) // conditional expression is constant
#endif

using namespace boost;

// See also test_conversion_from_string_to_integral(CharT)
// in libs/conversion/lexical_cast_test.cpp
template<class T, class CharT>
void test_too_long_number(CharT zero)
{
    typedef std::numeric_limits<T> limits;

    std::basic_string<CharT> s;

    std::basic_ostringstream<CharT> o;
    o << (limits::max)() << zero;
    s = o.str();
    BOOST_TEST_THROWS(lexical_cast<T>(s), bad_lexical_cast);
    s[s.size()-1] += static_cast<CharT>(9); // '0' -> '9'
    BOOST_TEST_THROWS(lexical_cast<T>(s), bad_lexical_cast);

    if (limits::is_signed)
    {
        std::basic_ostringstream<CharT> o2;
        o2 << (limits::min)() << zero;
        s = o2.str();
        BOOST_TEST_THROWS(lexical_cast<T>(s), bad_lexical_cast);
        s[s.size()-1] += static_cast<CharT>(9); // '0' -> '9'
        BOOST_TEST_THROWS(lexical_cast<T>(s), bad_lexical_cast);
    }
}

int main()
{
    test_too_long_number<boost::intmax_t>('0');
    test_too_long_number<boost::uintmax_t>('0');
#if !defined(BOOST_LCAST_NO_WCHAR_T)
    test_too_long_number<boost::intmax_t>(L'0');
    test_too_long_number<boost::uintmax_t>(L'0');
#endif

    return boost::report_errors();
}
