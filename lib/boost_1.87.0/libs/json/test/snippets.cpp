//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4101)
# pragma warning(disable: 4996)
#elif defined(__clang__)
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused"
#elif defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wunused"
#endif

#include <boost/json.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "test_suite.hpp"
#include "doc_types.hpp"

//[snippet_conv_spec_trait2
namespace boost
{
namespace json
{

template<>
struct is_sequence_like< user_ns::ip_address >
    : std::false_type
{ };

} // namespace json
} // namespace boost
//]

namespace user_ns2 {

class ip_address : public user_ns::ip_address
{
public:
    using user_ns::ip_address::ip_address;
};

using namespace boost::json;

//[snippet_tag_invoke_1
void
tag_invoke( const value_from_tag&, value& jv, ip_address const& addr )
{
    // Store the IP address as a 4-element array of octets
    const unsigned char* b = addr.begin();
    jv = { b[0], b[1], b[2], b[3] };
}

ip_address
tag_invoke( const value_to_tag< ip_address >&, value const& jv )
{
    array const& arr = jv.as_array();
    return ip_address(
        arr.at(0).to_number< unsigned char >(),
        arr.at(1).to_number< unsigned char >(),
        arr.at(2).to_number< unsigned char >(),
        arr.at(3).to_number< unsigned char >() );
}
//]

//[snippet_nothrow_1
result_for< ip_address, value >::type
tag_invoke( const try_value_to_tag< ip_address >&, value const& jv )
{
    if( !jv.is_array() )
        return make_error_code( std::errc::invalid_argument );

    array const& arr = jv.get_array();
    if( arr.size() != 4 )
        return make_error_code( std::errc::invalid_argument );

    boost::system::result< unsigned char > oct1
        = try_value_to< unsigned char >( arr[0] );
    if( !oct1 )
        return make_error_code( std::errc::invalid_argument );

    boost::system::result< unsigned char > oct2
        = try_value_to< unsigned char >( arr[1] );
    if( !oct2 )
        return make_error_code( std::errc::invalid_argument );

    boost::system::result< unsigned char > oct3
        = try_value_to< unsigned char >( arr[2] );
    if( !oct3 )
        return make_error_code( std::errc::invalid_argument );

    boost::system::result< unsigned char > oct4
        = try_value_to< unsigned char >( arr[3] );
    if( !oct4 )
        return make_error_code( std::errc::invalid_argument );

    return ip_address{ *oct1, *oct2, *oct3, *oct4 };
}
//]

} // namespace user_ns

//[doc_context_conversion_1
namespace user_ns
{

struct as_string
{ };

void
tag_invoke(
    boost::json::value_from_tag, boost::json::value& jv, const ip_address& addr, const as_string& )
{
    boost::json::string& js = jv.emplace_string();
    js.resize( 4 * 3 + 3 + 1 ); // XXX.XXX.XXX.XXX\0
    auto it = addr.begin();
    auto n = std::sprintf(
        js.data(), "%hhu.%hhu.%hhu.%hhu", it[0], it[1], it[2], it[3] );
    js.resize(n);
}

ip_address
tag_invoke(
    boost::json::value_to_tag<ip_address>, const boost::json::value& jv, const as_string& )
{
    const boost::json::string& js = jv.as_string();

    unsigned char octets[4];
    int result = std::sscanf(
        js.data(), "%hhu.%hhu.%hhu.%hhu", octets, octets + 1, octets + 2, octets + 3 );
    if( result != 4 )
        throw std::invalid_argument("not an IP address");

    return ip_address( octets[0], octets[1], octets[2], octets[3] );
}

}
//]

//[doc_context_conversion_4
namespace user_ns
{

struct as_iso_8601
{ };

void
tag_invoke(
    boost::json::value_from_tag, boost::json::value& jv, std::chrono::system_clock::time_point tp, const as_iso_8601& )
{
    boost::json::string& js = jv.emplace_string();
    js.resize( 4 + 2 * 5 + 5 + 1 ); // YYYY-mm-ddTHH:MM:ss\0

    std::time_t t = std::chrono::system_clock::to_time_t( tp );
    std::tm tm = *std::gmtime( &t );
    std::size_t n = std::strftime(
        js.data(), js.size(), "%FT%T", &tm );
    js.resize(n);
}

}
//]

//[doc_context_conversion_6
namespace user_ns
{

struct date_format
{
    std::string format;
    std::size_t buffer_size;
};

void
tag_invoke(
    boost::json::value_from_tag, boost::json::value& jv, std::chrono::system_clock::time_point tp, const date_format& ctx )
{
    boost::json::string& js = jv.emplace_string();
    js.resize( ctx.buffer_size );

    std::time_t t = std::chrono::system_clock::to_time_t( tp );
    std::size_t n = std::strftime(
        js.data(), js.size(), ctx.format.c_str(), std::gmtime( &t ) );
    js.resize(n);
}

}
//]

//[doc_context_conversion_10
namespace user_ns
{

struct maps_as_objects
{ };

template<
    class Key,
    class Value,
    class Ctx >
void
tag_invoke(
    boost::json::value_from_tag,
    boost::json::value& jv,
    const std::map<Key, Value>& m,
    const maps_as_objects&,
    const Ctx& ctx )
{
    boost::json::object& jo = jv.emplace_object();

    for( const auto& item: m )
    {
        auto k = boost::json::value_from( item.first, ctx, jo.storage() );
        auto v = boost::json::value_from( item.second, ctx, jo.storage() );
        jo[std::move( k.as_string() )] = std::move( v );
    }
}

template<
    class Key,
    class Value,
    class Ctx >
std::map<Key, Value>
tag_invoke(
    boost::json::value_to_tag< std::map<Key, Value> >,
    boost::json::value const& jv,
    const maps_as_objects&,
    const Ctx& ctx )
{
    const boost::json::object& jo = jv.as_object();
    std::map< Key, Value > result;
    for( auto&& item: jo )
    {
        Key k = boost::json::value_to< Key >( item.key(), ctx );
        Value v = boost::json::value_to< Value >( item.value(), ctx );
        result.emplace( std::move(k), std::move(v) );
    }
    return result;
}

}
//]

namespace boost {
namespace json {

namespace {

//[snippet_strings_5
string greeting( string_view first_name, string_view last_name )
{
    const char hello[] = "Hello, ";
    const std::size_t sz = first_name.size() + last_name.size() + sizeof(hello) + 1;

    string js;
    js.reserve(sz);

    char* p = std::copy( hello, hello + sizeof(hello) - 1, js.data() );
    p = std::copy( first_name.begin(), first_name.end(), p );
    *p++ = ' ';
    p = std::copy( last_name.begin(), last_name.end(), p );
    *p++ = '!';

    js.grow( sz );
    return js;
}
//]

void
usingStrings()
{
    {
        //[snippet_strings_1

        string str1; // empty string, uses the default memory resource

        string str2( make_shared_resource<monotonic_resource>() ); // empty string, uses a counted monotonic resource

        //]
    }
    {
        //[snippet_strings_2

        std::string sstr1 = "helloworld";
        std::string sstr2 = "world";

        json::string jstr1 = "helloworld";
        json::string jstr2 = "world";

        assert( jstr2.insert(0, jstr1.subview(0, 5)) == "helloworld" );

        // this is equivalent to
        assert( sstr2.insert(0, sstr1, 0, 5) == "helloworld" );

        //]
    }
    {
        //[snippet_strings_3

        std::string sstr = "hello";

        json::string jstr = "hello";

        assert(sstr.append({'w', 'o', 'r', 'l', 'd'}) == "helloworld");

        // such syntax is inefficient, and the same can
        // be achieved with a character array.

        assert(jstr.append("world") == "helloworld");

        //]
    }

    {
        //[snippet_strings_4

        json::string str = "Boost.JSON";
        json::string_view sv = str;

        // all of these call compare(string_view)
        str.compare(sv);

        str.compare(sv.substr(0, 5));

        str.compare(str);

        str.compare("Boost");

        //]
    }

    {
        auto str = greeting( "John", "Smith" );
        (void)str;
        assert( str == "Hello, John Smith!");
    }
}

//----------------------------------------------------------

void
usingValues()
{
    {
        //[snippet_value_1

        value jv1;
        value jv2( nullptr );

        assert( jv1.is_null() );
        assert( jv2.is_null() );

        //]
    }
    {
        //[snippet_value_2

        value jv( object_kind );

        assert( jv.kind() == kind::object );
        assert( jv.is_object() );
        assert( ! jv.is_number() );

        //]
    }
    {
        auto f = []{
        //[snippet_value_3

        value jv( object_kind );

        if( auto p = jv.if_object() )
            return p->size();

        //]
        return std::size_t(0);
        };
        (void)f;
    }
    {
        //[snippet_value_4

        value jv;
        jv = value( array_kind );

        assert( jv.is_array() );

        jv.emplace_string();

        assert( jv.is_string() );

        //]
    }
    {
        //[snippet_value_5

        value jv;
        jv.emplace_string() = "Hello, world!";

        int64_t& num = jv.emplace_int64();
        num = 1;

        assert( jv.is_int64() );

        //]
    }
    {
        try
        {
            //[snippet_value_6

            value jv( true );
            jv.as_bool() = true;

            jv.as_string() = "Hello, world!"; // throws an exception

            //]
        }
        catch(...)
        {
        }
    }
    {
        //[snippet_value_7

        value jv( string_kind );
        if( string* str = jv.if_string() )
            *str = "Hello, world!";

        //]
    }
    {
        //[snippet_value_8

        value jv( string_kind );

        // The compiler's static analysis can see that
        // a null pointer is never dereferenced.
        *jv.if_string() = "Hello, world!";

        //]
    }
    {
        //[snippet_value_9

        value jv( string_kind );
        if( boost::system::result<string&> str = jv.try_as_string() )
            *str = "Hello, world!";

        try
        {
            jv.try_as_bool().value() = true;
        }
        catch(...)
        {
        }

        //]
    }
}

//----------------------------------------------------------

void
usingInitLists()
{
    {
        //[snippet_init_list_1

        value jv = {
            { "name", "John Doe" },
            { "active", true },
            { "associated-accounts", nullptr },
            { "total-balance", 330.00 },
            { "account-balances", { 84, 120, 126 } } };

        //]
    }

    {
        //[snippet_init_list_2

        value jv = { true, 2, "hello", nullptr };

        assert( jv.is_array() );

        assert( jv.as_array().size() == 4 );

        assert( serialize(jv) == R"([true,2,"hello",null])" );

        //]
    }

    {
        //[snippet_init_list_3

        value jv = { true, 2, "hello", { "bye", nullptr, false } };

        assert( jv.is_array() );

        assert( jv.as_array().back().is_array() );

        assert( serialize(jv) == R"([true,2,"hello",["bye",null,false]])" );

        //]
    }

    {
        //[snippet_init_list_4

        // Should this be an array or an object?
        value jv = { { "hello", 42 }, { "world", 43 } };

        //]
    }

    {
        //[snippet_init_list_5

        value jv1 = { { "hello", 42 }, { "world", 43 } };

        assert( jv1.is_object() );

        assert( jv1.as_object().size() == 2 );

        assert( serialize(jv1) == R"({"hello":42,"world":43})" );

        // All of the following are arrays

        value jv2 = { { "make", "Tesla" }, { "model", 3 }, "black" };

        value jv3 = { { "library", "JSON" }, { "Boost", "C++", "Fast", "JSON" } };

        value jv4 = { { "color", "blue" }, { 1, "red" } };

        assert( jv2.is_array() && jv3.is_array() && jv4.is_array() );

        //]
    }

    {
        //[snippet_init_list_6

        value jv = { { "hello", 42 }, array{ "world", 43 } };

        assert( jv.is_array() );

        array& ja = jv.as_array();

        assert( ja[0].is_array() && ja[1].is_array());

        assert ( serialize(jv) == R"([["hello",42],["world",43]])" );

        //]

        (void)ja;
    }

    {
        //[snippet_init_list_7

        value jv = { { "mercury", 36 }, { "venus", 67 }, { "earth", 93 } };

        assert( jv.is_object() );

        assert( serialize(jv) == R"({"mercury":36,"venus":67,"earth":93})" );

        array ja = { { "mercury", 36 }, { "venus", 67 }, { "earth", 93 } };

        assert( serialize(ja) == R"([["mercury",36],["venus",67],["earth",93]])" );

        //]

        (void)ja;
    }

    {
        //[snippet_init_list_8

        object jo = { { "mercury", { { "distance", 36 } } }, { "venus", { 67, "million miles" } }, { "earth", 93 } };

        assert( jo["mercury"].is_object() );

        assert( jo["venus"].is_array() );

        //]
    }

    {
        //[snippet_init_list_9

        object jo1 = { { "john", 100 }, { "dave", 500 }, { "joe", 300 } };

        value jv = { { "clients", std::move(jo1) } };

        object& jo2 = jv.as_object()["clients"].as_object();

        assert( ! jo2.empty() && jo1.empty() );

        assert( serialize(jv) == R"({"clients":{"john":100,"dave":500,"joe":300}})" );

        //]

        (void)jo2;
    }
}

//----------------------------------------------------------

void
usingArrays()
{
    {
        //[snippet_arrays_1

        array arr1; // empty array, uses the default memory resource

        array arr2( make_shared_resource<monotonic_resource>() ); // empty array, uses a counted monotonic resource

        //]
    }
    {
        //[snippet_arrays_2

        array arr( { "Hello", 42, true } );

        //]
    }
    try
    {
        //[snippet_arrays_3

        array arr;

        arr.emplace_back( "Hello" );
        arr.emplace_back( 42 );
        arr.emplace_back( true );

        //]

        //[snippet_arrays_4

        assert( arr[0].as_string() == "Hello" );

        // The following line throws system_error, since the index is out of range
        arr.at( 3 ) = nullptr;

        //]
    }
    catch (...)
    {
    }
}

//----------------------------------------------------------

void
usingObjects()
{
    {
        //[snippet_objects_1

        object obj1; // empty object, uses the default memory resource

        object obj2( make_shared_resource<monotonic_resource>() ); // empty object, uses a counted monotonic resource

        //]
    }
    {
        //[snippet_objects_2

        object obj( {{"key1", "value1" }, { "key2", 42 }, { "key3", false }} );

        //]
    }
    {
        //[snippet_objects_3

        object obj;

        obj.emplace( "key1", "value1" );
        obj.emplace( "key2", 42 );
        obj.emplace( "key3", false );

        //]
    }
    try
    {
        //[snippet_objects_4

        object obj;

        obj["key1"] = "value1";
        obj["key2"] = 42;
        obj["key3"] = false;

        // The following line throws system_error, since the key does not exist
        obj.at( "key4" );

        //]
    }
    catch (...)
    {
    }
    {
        //[snippet_objects_5

        object obj{{"arr", {1, 11}}};
        value& arr = obj.at("arr");
        obj.emplace("added", "value"); // invalidates arr

        //]

        (void)arr;
    }
}

//[snippet_conv_5

template< class T >
struct vec3
{
    T x, y, z;
};

template< class T >
void tag_invoke( const value_from_tag&, value& jv, const vec3<T>& vec )
{
    jv = {
        { "x", vec.x },
        { "y", vec.y },
        { "z", vec.z }
    };
}

//]

//[snippet_conv_10

struct customer
{
    std::uint64_t id;
    std::string name;
    bool late;

    customer() = default;

    customer( std::uint64_t i, const std::string& n, bool l )
        : id( i ), name( n ), late( l ) { }
};

void tag_invoke( const value_from_tag&, value& jv, customer const& c )
{
    // Assign a JSON value
    jv = {
        { "id", c.id },
        { "name", c.name },
        { "late", c.late }
    };
}

//]

//[snippet_conv_14

customer tag_invoke( const value_to_tag<customer>&, const value& jv )
{
    // at() throws if `jv` is not an object, or if the key is not found.
    // as_uint64() will throw if the value is not an unsigned 64-bit integer.
    std::uint64_t id = jv.at( "id" ).as_uint64();

    // We already know that jv is an object from
    // the previous call to jv.as_object() succeeding,
    // now we use jv.get_object() which skips the
    // check. value_to will throw if jv.kind() != kind::string
    std::string name = value_to< std::string >( jv.get_object().at( "name" ) );

    // id and name are constructed from JSON in the member
    // initializer list above, but we can also use regular
    // assignments in the body of the function as shown below.
    // as_bool() will throw if kv.kind() != kind::bool
    bool late = jv.get_object().at( "late" ).as_bool();

    return customer(id, name, late);
}

//]

void
usingExchange()
{
    {
        //[snippet_conv_1

        std::vector< int > v1{ 1, 2, 3, 4 };

        // Convert the vector to a JSON array
        value jv = value_from( v1 );
        assert( serialize( jv ) == R"([1,2,3,4])" );

        // Convert back to vector< int >
        std::vector< int > v2 = value_to< std::vector< int > >( jv );
        assert( v1 == v2 );

        //]

        (void)v2;
    }
    {
        using namespace user_ns2;

        //[snippet_tag_invoke_3
        std::map< std::string, ip_address > computers = {
            { "Alex", { 192, 168, 1, 1 } },
            { "Blake", { 192, 168, 1, 2 } },
            { "Carol", { 192, 168, 1, 3 } },
        };

        // conversions are applied recursively;
        // the key type and value type will be converted
        // using value_from as well
        value jv = value_from( computers );
        assert( jv.is_object() );

        value serialized = parse(R"(
            {
                "Alex":  [ 192, 168, 1, 1 ],
                "Blake": [ 192, 168, 1, 2 ],
                "Carol": [ 192, 168, 1, 3 ]
            }
            )");
        assert( jv == serialized );
        //]

        (void)jv;
    }
    {
        using namespace user_ns2;

        //[snippet_tag_invoke_2
        ip_address addr = { 127, 0, 0, 12 };
        value jv = value_from( addr );
        assert( serialize( jv ) == R"([127,0,0,12])" );

        // Convert back to IP address
        ip_address addr2 = value_to< ip_address >( jv );
        assert(std::equal(
            addr.begin(), addr.end(), addr2.begin() ));
        //]

        (void)addr2;
    }
    {
        using namespace user_ns2;

        //[snippet_nothrow_2
        value jv = parse( R"([127,0,0,12])" );

        boost::system::result< ip_address > addr
            = try_value_to< ip_address >( jv );
        assert( addr.has_value() );

        ip_address addr2{ 127, 0, 0, 12 };
        assert(std::equal(
            addr->begin(), addr->end(), addr2.begin() ));

        // this fails without exception
        addr = try_value_to< ip_address >( value() );
        assert( addr.has_error() );
        //]

        (void)addr;
        (void)addr2;
    }
    {
        //[snippet_conv_recursive
        std::map< std::string, std::pair<int, bool> > m = {
            {"a", {1, false}},
            {"b", {4, true}},
            {"c", {5, false}},
        };

        value jv = value_from( m );

        assert(( jv == object{
            {"a", array{1, false}},
            {"b", array{4, true}},
            {"c", array{5, false}},
        }));
        //]
    }
}

void
usingPointer()
{
    //[snippet_pointer_1
    value jv = { {"one", 1}, {"two", 2} };
    assert( jv.at("one") == jv.at_pointer("/one") );

    jv.at_pointer("/one") = {{"foo", "bar"}};
    assert( jv.at("one").at("foo") == jv.at_pointer("/one/foo") );

    jv.at_pointer("/one/foo") = {true, 4, "qwerty"};
    assert( jv.at("one").at("foo").at(1) == jv.at_pointer("/one/foo/1") );
    //]

    value* elem1 = [&]() -> value*
    {
        //[snippet_pointer_2
        object* obj = jv.if_object();
        if( !obj )
            return nullptr;

        value* val = obj->if_contains("one");
        if( !val )
            return nullptr;

        obj = val->if_object();
        if( !obj )
            return nullptr;

        val = obj->if_contains("foo");
        if( !val )
            return nullptr;

        array* arr = val->if_array();
        if( !arr )
            return nullptr;

        return arr->if_contains(1);
        //]
    }();

    value* elem2 = [&]() -> value*
    {
        //[snippet_pointer_3
        boost::system::error_code ec;
        return jv.find_pointer("/one/foo/1", ec);
        //]
    }();

    (void)elem1;
    (void)elem2;
    assert( elem1 == elem2 );
}


void
usingSetAtPointer()
{
    //[snippet_pointer_4
    value jv;
    jv.set_at_pointer("/two/bar/0", 12);
    assert( jv.is_object() );
    assert( jv.at_pointer("/two").is_object() );
    assert( jv.at_pointer("/two/bar").is_array() );
    assert( jv.at_pointer("/two/bar/0") == 12 );
    //]

    jv = nullptr;
    //[snippet_pointer_5
    set_pointer_options opts;
    opts.create_arrays = false;

    jv.set_at_pointer("/two/bar/0", 12, opts);
    assert( jv.is_object() );
    assert( jv.at_pointer("/two").is_object() );
    assert( jv.at_pointer("/two/bar").is_object() ); // object, not array
    assert( jv.at_pointer("/two/bar/0") == 12 );
    //]
}

BOOST_STATIC_ASSERT(
    has_value_from<customer>::value);

BOOST_STATIC_ASSERT(
    has_value_from<user_ns2::ip_address>::value);
BOOST_STATIC_ASSERT(
    has_value_to<user_ns2::ip_address>::value);

} // (anon)

} // json
} // boost

//----------------------------------------------------------


namespace {

class my_non_deallocating_resource { };

} // (anon)

//[snippet_allocators_14
namespace boost {
namespace json {

template<>
struct is_deallocate_trivial< my_non_deallocating_resource >
{
    static constexpr bool value = true;
};

} // json
} // boost

//]
namespace boost {
namespace json {

namespace
{

void
usingSpecializedTrait()
{
    value jv1{127, 0, 0, 1};
    auto const addr = value_to< user_ns::ip_address >( jv1 );
    auto const jv2 = value_from(addr);
    assert( jv1 == jv2 );
}

void
usingContextualConversions()
{
    using namespace user_ns;
    using namespace boost::json;
    {
//[doc_context_conversion_2
        ip_address addr( 192, 168, 10, 11 );

        value jv = value_from( addr, as_string() );
        assert( jv == parse(R"( "192.168.10.11" )") );

        ip_address addr2 = value_to< ip_address >( jv, as_string() );
        assert(std::equal(
            addr.begin(), addr.end(), addr2.begin() ));
//]
        (void)addr2;
    }

    {
//[doc_context_conversion_3
        std::map< std::string, ip_address > computers = {
            { "Alex", { 192, 168, 1, 1 } },
            { "Blake", { 192, 168, 1, 2 } },
            { "Carol", { 192, 168, 1, 3 } },
        };
        value jv = value_from( computers, as_string() );
        assert( jv == parse(
            "{                               "
            "    \"Alex\" : \"192.168.1.1\", "
            "    \"Blake\": \"192.168.1.2\", "
            "    \"Carol\": \"192.168.1.3\"  "
            "}                               "
            ) );
//]
        (void)jv;
    }

    {
//[doc_context_conversion_5
        std::chrono::system_clock::time_point tp;
        value jv = value_from( tp, as_iso_8601() );
        assert( jv == parse(R"( "1970-01-01T00:00:00" )") );
//]
        (void)jv;
    }

    {
//[doc_context_conversion_7
        std::chrono::system_clock::time_point tp;

        value jv = value_from( tp, date_format{ "%T %D", 18 } );
        assert( jv == parse(R"( "00:00:00 01/01/70" )") );

        jv = value_from( tp, as_iso_8601() );
        assert( jv == parse(R"( "1970-01-01T00:00:00" )") );
//]
        (void)jv;
    }

    {
//[doc_context_conversion_8
        using time_point = std::chrono::system_clock::time_point;
        time_point start;
        std::vector< std::pair<time_point, ip_address> > log = {
            { start += std::chrono::seconds(10), {192, 168, 10, 11} },
            { start += std::chrono::hours(2),    {192, 168, 10, 13} },
            { start += std::chrono::minutes(14), {192, 168, 10, 10} },
        };
        value jv = value_from(
            log, std::make_tuple( as_string(), as_iso_8601() ) );
        assert( jv == parse(
            " [                                                   "
            "     [ \"1970-01-01T00:00:10\", \"192.168.10.11\" ], "
            "     [ \"1970-01-01T02:00:10\", \"192.168.10.13\" ], "
            "     [ \"1970-01-01T02:14:10\", \"192.168.10.10\" ]  "
            " ]                                                   "
            ) );
//]
        (void)jv;
    }

    {
        using time_point = std::chrono::system_clock::time_point;
        time_point start;
//[doc_context_conversion_9

        std::map< time_point, ip_address > log = {
            { start += std::chrono::seconds(10), {192, 168, 10, 11} },
            { start += std::chrono::hours(2),    {192, 168, 10, 13} },
            { start += std::chrono::minutes(14), {192, 168, 10, 10} },
        };

        value jv = value_from(
            log,
            std::make_tuple( maps_as_objects(), as_string(), as_iso_8601() ) );
        assert( jv == parse(
            " {                                               "
            "     \"1970-01-01T00:00:10\": \"192.168.10.11\", "
            "     \"1970-01-01T02:00:10\": \"192.168.10.13\", "
            "     \"1970-01-01T02:14:10\": \"192.168.10.10\"  "
            " }                                               "
            ) );
//]
        (void)jv;
    }
}

void
usingParseInto()
{
//[doc_parse_into_1
    std::map< std::string, std::vector<int> > vectors;
    string_view input = R"( { "even": [2,4,6], "odd": [1,3,5] } )";
    parse_into(vectors, input);
//]

    std::string output = serialize(vectors);
    (void)output;
    assert( output == R"({"even":[2,4,6],"odd":[1,3,5]})" );
}

} // namespace

class snippets_test
{
public:
    void
    run()
    {
        usingValues();
        usingInitLists();
        usingExchange();
        usingArrays();
        usingObjects();
        usingStrings();
        usingPointer();
        usingSpecializedTrait();
        usingSetAtPointer();
        usingContextualConversions();
        usingParseInto();

        BOOST_TEST_PASS();
    }
};

TEST_SUITE(snippets_test, "boost.json.snippets");

} // namespace json
} // namespace boost
