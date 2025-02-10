// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config.hpp>

#if defined(BOOST_GCC) && BOOST_GCC < 50000
# define BOOST_ALLOW_DEPRECATED
#endif

#include <boost/variant2/variant.hpp>
#include <boost/json/value_from.hpp>
#include <boost/json/serialize.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>

using namespace boost::variant2;
namespace json = boost::json;

int main()
{
    {
        monostate m;
        json::value w = json::value_from( m );
        BOOST_TEST_EQ( w, json::value( nullptr ) );
    }

    {
        variant<monostate, int, std::string> v;
        json::value w = json::value_from( v );
        BOOST_TEST_EQ( w, json::value( nullptr ) );
    }

    {
        variant<monostate, int, std::string> v( 17 );
        json::value w = json::value_from( v );
        BOOST_TEST_EQ( w, json::value( 17 ) );
    }

    {
        variant<monostate, int, std::string> v( "test" );
        json::value w = json::value_from( v );
        BOOST_TEST_EQ( w, json::value( "test" ) );
    }

    return boost::report_errors();
}
