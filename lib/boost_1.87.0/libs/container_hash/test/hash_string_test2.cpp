// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/utility/string_view.hpp>
#include <boost/core/detail/string_view.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <string>
#if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
#include <string_view>
#endif
#include <vector>
#include <deque>
#include <list>
#if !defined(BOOST_NO_CXX11_HDR_FORWARD_LIST)
# include <forward_list>
#endif

// Test whether the hash values of a string and a
// string_view that refers to it match. This is
// important for unordered heterogeneous lookups.

template<class T> std::size_t hv( T const& t )
{
    return boost::hash<T>()( t );
}

int main()
{
    std::string s( "Test." );
    std::size_t h0 = hv( s );

    BOOST_TEST_EQ( h0, hv( boost::string_view( s ) ) );
    BOOST_TEST_EQ( h0, hv( boost::core::string_view( s ) ) );

#if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
    BOOST_TEST_EQ( h0, hv( std::string_view( s ) ) );
#endif

    BOOST_TEST_EQ( h0, hv( std::vector<char>( s.begin(), s.end() ) ) );
    BOOST_TEST_EQ( h0, hv( std::deque<char>( s.begin(), s.end() ) ) );
    BOOST_TEST_EQ( h0, hv( std::list<char>( s.begin(), s.end() ) ) );

#if !defined(BOOST_NO_CXX11_HDR_FORWARD_LIST)
    BOOST_TEST_EQ( h0, hv( std::forward_list<char>( s.begin(), s.end() ) ) );
#endif

    return boost::report_errors();
}
