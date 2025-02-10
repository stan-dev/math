//  Copyright Antony Polukhin, 2012-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <boost/lexical_cast/detail/converter_lexical.hpp>

#include <boost/core/lightweight_test.hpp>

#include <boost/array.hpp>
#include <boost/container/string.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/utility/string_view.hpp>

template <class T>
struct is_optimized_stream : boost::false_type {};

template <class CharT, class Traits, std::size_t CharacterBufferSize>
struct is_optimized_stream< boost::detail::lcast::optimized_src_stream<CharT, Traits, CharacterBufferSize> > : boost::true_type {};

template <class Source, class Target>
static void assert_optimized_stream()
{
    BOOST_TEST((is_optimized_stream<
        typename boost::detail::lexical_converter_impl<Source, Target>::from_src_stream
    >::value));
}

template <class T>
static void test_optimized_types_to_string_const()
{
    namespace de = boost::detail;
    typedef de::lexical_cast_stream_traits<T, std::string> trait_1;
    BOOST_TEST((boost::is_same<typename trait_1::src_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_1::target_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_1::char_type, char>::value));
    assert_optimized_stream<T, std::string>();
    assert_optimized_stream<T, boost::container::string>();

    typedef de::lexical_cast_stream_traits<const T, std::string> trait_2;
    BOOST_TEST((boost::is_same<typename trait_2::src_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_2::target_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_2::char_type, char>::value));

    typedef de::lexical_cast_stream_traits<T, std::wstring> trait_3;
    BOOST_TEST((boost::is_same<typename trait_3::src_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_3::target_char_t, wchar_t>::value));
    BOOST_TEST((boost::is_same<typename trait_3::char_type, wchar_t>::value));
    assert_optimized_stream<T, std::wstring>();
    assert_optimized_stream<T, boost::container::wstring>();
}


template <class T>
static void test_optimized_types_to_string()
{
    test_optimized_types_to_string_const<T>();

    namespace de = boost::detail;
    typedef de::lexical_cast_stream_traits<std::string, T> trait_4;
    BOOST_TEST((boost::is_same<typename trait_4::src_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_4::target_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_4::char_type, char>::value));
    assert_optimized_stream<std::string, T>();

    typedef de::lexical_cast_stream_traits<const std::string, T> trait_5;
    BOOST_TEST((boost::is_same<typename trait_5::src_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_5::target_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_5::char_type, char>::value));

    typedef de::lexical_cast_stream_traits<const std::wstring, T> trait_6;
    BOOST_TEST((boost::is_same<typename trait_6::src_char_t, wchar_t>::value));
    BOOST_TEST((boost::is_same<typename trait_6::target_char_t, char>::value));
    BOOST_TEST((boost::is_same<typename trait_6::char_type, wchar_t>::value));
}


template <class T>
static void test_optimized_types_to_x_string()
{
    test_optimized_types_to_string<T>();
    assert_optimized_stream<std::wstring, T>();
    assert_optimized_stream<boost::container::wstring, T>();
}


template <class T>
static void test_optimized_types_to_wstring()
{
    assert_optimized_stream<std::wstring, T>();
    assert_optimized_stream<T, std::wstring>();
    assert_optimized_stream<boost::container::wstring, T>();
    assert_optimized_stream<T, boost::container::wstring>();
}

void test_metafunctions()
{
    test_optimized_types_to_x_string<bool>();
    test_optimized_types_to_x_string<char>();
    test_optimized_types_to_x_string<unsigned char>();
    test_optimized_types_to_x_string<signed char>();
    test_optimized_types_to_x_string<short>();
    test_optimized_types_to_x_string<unsigned short>();
    test_optimized_types_to_x_string<int>();
    test_optimized_types_to_x_string<unsigned int>();
    test_optimized_types_to_x_string<long>();
    test_optimized_types_to_x_string<unsigned long>();

#if defined(BOOST_HAS_LONG_LONG)
    test_optimized_types_to_x_string<boost::ulong_long_type>();
    test_optimized_types_to_x_string<boost::long_long_type>();
#elif defined(BOOST_HAS_MS_INT64)
    test_optimized_types_to_x_string<unsigned __int64>();
    test_optimized_types_to_x_string<__int64>();
#endif

    test_optimized_types_to_string<float>();
    test_optimized_types_to_string<double>();
    test_optimized_types_to_string<long double>();
    test_optimized_types_to_string<std::string>();
    test_optimized_types_to_string<char*>();
    //test_optimized_types_to_string<char[5]>();
    //test_optimized_types_to_string<char[1]>();
    test_optimized_types_to_string<unsigned char*>();
    //test_optimized_types_to_string<unsigned char[5]>();
    //test_optimized_types_to_string<unsigned char[1]>();
    test_optimized_types_to_string<signed char*>();
    //test_optimized_types_to_string<signed char[5]>();
    //test_optimized_types_to_string<signed char[1]>();
    test_optimized_types_to_string<boost::array<char, 1> >();
    test_optimized_types_to_string<boost::array<char, 5> >();
    test_optimized_types_to_string<boost::array<unsigned char, 1> >();
    test_optimized_types_to_string<boost::array<unsigned char, 5> >();
    test_optimized_types_to_string<boost::array<signed char, 1> >();
    test_optimized_types_to_string<boost::array<signed char, 5> >();
    test_optimized_types_to_string<boost::iterator_range<char*> >();
    test_optimized_types_to_string<boost::iterator_range<unsigned char*> >();
    test_optimized_types_to_string<boost::iterator_range<signed char*> >();

    test_optimized_types_to_string_const<boost::array<const char, 1> >();
    test_optimized_types_to_string_const<boost::array<const char, 5> >();
    test_optimized_types_to_string_const<boost::array<const unsigned char, 1> >();
    test_optimized_types_to_string_const<boost::array<const unsigned char, 5> >();
    test_optimized_types_to_string_const<boost::array<const signed char, 1> >();
    test_optimized_types_to_string_const<boost::array<const signed char, 5> >();
    test_optimized_types_to_string_const<boost::iterator_range<const char*> >();
    test_optimized_types_to_string_const<boost::iterator_range<const unsigned char*> >();
    test_optimized_types_to_string_const<boost::iterator_range<const signed char*> >();

#ifndef BOOST_NO_CXX11_HDR_ARRAY
    test_optimized_types_to_string<std::array<char, 1> >();
    test_optimized_types_to_string<std::array<char, 5> >();
    test_optimized_types_to_string<std::array<unsigned char, 1> >();
    test_optimized_types_to_string<std::array<unsigned char, 5> >();
    test_optimized_types_to_string<std::array<signed char, 1> >();
    test_optimized_types_to_string<std::array<signed char, 5> >();

    test_optimized_types_to_string_const<std::array<const char, 1> >();
    test_optimized_types_to_string_const<std::array<const char, 5> >();
    test_optimized_types_to_string_const<std::array<const unsigned char, 1> >();
    test_optimized_types_to_string_const<std::array<const unsigned char, 5> >();
    test_optimized_types_to_string_const<std::array<const signed char, 1> >();
    test_optimized_types_to_string_const<std::array<const signed char, 5> >();

    test_optimized_types_to_wstring<std::array<wchar_t, 42>>();
    test_optimized_types_to_wstring<std::array<const wchar_t, 42>>();
#endif
    
    test_optimized_types_to_wstring<wchar_t*>();
    test_optimized_types_to_wstring<const wchar_t*>();
    test_optimized_types_to_wstring<boost::iterator_range<const wchar_t*>>();
    test_optimized_types_to_wstring<boost::iterator_range<wchar_t*>>();

    test_optimized_types_to_string<boost::string_view>();
    test_optimized_types_to_wstring<boost::basic_string_view<wchar_t>>();

#ifndef BOOST_NO_CXX17_HDR_STRING_VIEW
    test_optimized_types_to_string<std::string_view>();
    test_optimized_types_to_wstring<std::basic_string_view<wchar_t>>();
#endif
}

int main()
{
    test_metafunctions();

    return boost::report_errors();
}

