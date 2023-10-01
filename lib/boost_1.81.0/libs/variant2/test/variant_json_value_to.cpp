// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/variant2/variant.hpp>
#include <boost/json/value_to.hpp>
#include <boost/json/serialize.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::variant2;
namespace json = boost::json;

int main()
{
    {
        json::value v;
        auto r = json::try_value_to<monostate>( v );
        BOOST_TEST( r.has_value() );
    }

    using V = variant<monostate, int, std::string>;

    {
        json::value v;
        auto r = json::try_value_to<V>( v );
        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( *r, V() );
    }

    {
        json::value v( 12 );
        auto r = json::try_value_to<V>( v );
        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( *r, V(12) );
    }

    {
        json::value v( "test" );
        auto r = json::try_value_to<V>( v );
        BOOST_TEST( r.has_value() ) && BOOST_TEST_EQ( *r, V("test") );
    }

    return boost::report_errors();
}
