//
// Copyright (c) 2021 Dmitry Arkhipov (grisumbras@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

//[doc_forward_conversion_1
namespace boost {
namespace json {

class value;

struct value_from_tag;

template< class T >
struct try_value_to_tag;

template< class T1, class T2 >
struct result_for;

template< class T >
void value_from( T&& t, value& jv );

template< class T >
typename result_for< T, value >::type
try_value_to( const value& jv );

}
}
//]

#include "doc_types.hpp"

#include <system_error>

//[doc_forward_conversion_2
namespace user_ns
{

template< class JsonValue >
void tag_invoke( const boost::json::value_from_tag&, JsonValue& jv, const ip_address& addr )
{
    const unsigned char* b = addr.begin();
    jv = { b[0], b[1], b[2], b[3] };
}

template< class JsonValue >
typename boost::json::result_for< ip_address, JsonValue >::type
tag_invoke(
    const boost::json::try_value_to_tag< ip_address >&,
    const JsonValue& jv )
{
    using namespace boost::json;

    if( !jv.is_array() )
        return make_error_code( std::errc::invalid_argument );

    auto const& arr = jv.get_array();
    if( arr.size() != 4 )
        return make_error_code( std::errc::invalid_argument );

    auto oct1 = try_value_to< unsigned char >( arr[0] );
    if( !oct1 )
        return make_error_code( std::errc::invalid_argument );

    auto oct2 = try_value_to< unsigned char >( arr[1] );
    if( !oct2 )
        return make_error_code( std::errc::invalid_argument );

    auto oct3 = try_value_to< unsigned char >( arr[2] );
    if( !oct3 )
        return make_error_code( std::errc::invalid_argument );

    auto oct4 = try_value_to< unsigned char >( arr[3] );
    if( !oct4 )
        return make_error_code( std::errc::invalid_argument );

    return ip_address{ *oct1, *oct2, *oct3, *oct4 };
}

}
//]


#include <boost/json/value_from.hpp>
#include <boost/json/value_to.hpp>

#include "test_suite.hpp"


BOOST_JSON_NS_BEGIN

class doc_forward_conversion
{
public:
    void
    run()
    {
        value const jv{ 212, 115, 81, 22 };
        auto const addr = value_to< user_ns::ip_address >( jv );
        BOOST_TEST( get<0>(addr) == 212 );
        BOOST_TEST( get<1>(addr) == 115 );
        BOOST_TEST( get<2>(addr) == 81 );
        BOOST_TEST( get<3>(addr) == 22 );

        value const jv2 = value_from( addr );
        BOOST_TEST( jv == jv2 );
    }
};

TEST_SUITE(doc_forward_conversion, "boost.json.doc_forward_conversion");

BOOST_JSON_NS_END
