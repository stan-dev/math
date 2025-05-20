//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2021 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/parse_into.hpp>

#include <boost/json/serialize.hpp>
#include <boost/json/value_from.hpp>
#include <boost/json/value_to.hpp>
#include <boost/describe.hpp>

#include <climits>
#include <map>

#ifndef BOOST_NO_CXX17_HDR_FILESYSTEM
#include <filesystem>
#endif // BOOST_NO_CXX17_HDR_FILESYSTEM

#include "test.hpp"
#include "test_suite.hpp"

struct X
{
    int a;
    float b;
    std::string c;
};

BOOST_DESCRIBE_STRUCT(X, (), (a, b, c))

bool operator==( X const& x1, X const& x2 )
{
    return x1.a == x2.a && x1.b == x2.b && x1.c == x2.c;
}

struct Y
{
    std::vector<X> v;
    std::map<std::string, X> m;
};

BOOST_DESCRIBE_STRUCT(Y, (), (v, m))

bool operator==( Y const& y1, Y const& y2 )
{
    return y1.v == y2.v && y1.m == y2.m;
}

struct Z : X
{
    Z() = default;
    Z(X x, bool b) : X(x), d(b)
    {}

private:
    bool d = false;

    friend
    bool operator==( Z const& z1, Z const& z2 )
    {
        X const& x1 = z1;
        X const& x2 = z2;
        return (x1 == x2) && z1.d == z2.d;
    }

    BOOST_DESCRIBE_CLASS(Z, (X), (), (), (d))
};

BOOST_DEFINE_ENUM_CLASS(E, x, y, z)

namespace boost {
namespace json {

template<>
struct is_described_class<Z>
    : std::true_type
{ };

class parse_into_test
{
public:

    template<class T>
    void testParseIntoValue( value const& jv )
    {
#if defined(__GNUC__) && __GNUC__ < 5
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
        T t1 = value_to<T>(jv);
        (void)t1; // older GCC thinks t1 can be unused
        std::string json = serialize(jv);

        T t2{};
        system::error_code jec;
        parse_into(t2, json, jec);
        BOOST_TEST( !jec.failed() ) && BOOST_TEST( t1 == t2 );

        T t3{};
        std::error_code ec;
        parse_into(t3, json, ec);
        BOOST_TEST( !ec ) && BOOST_TEST( t1 == t3 );

        T t4{};
        parse_into(t4, json);
        BOOST_TEST( t1 == t4 );

        std::istringstream is(json);
        T t5{};
        jec = {};
        parse_into(t5, is, jec);
        BOOST_TEST( !jec.failed() ) && BOOST_TEST( t1 == t5 );

        is.clear();
        is.seekg(0);
        T t6{};
        ec = {};
        parse_into(t6, is, ec);
        BOOST_TEST( !ec ) && BOOST_TEST( t1 == t6 );

        is.str(json);
        is.clear();
        T t7{};
        parse_into(t7, is);
        BOOST_TEST( t1 == t7 );

        parse_options opt;
        opt.allow_comments = true;
        json = "// this is a comment\n" + json;

        T t8{};
        parser_for<T> p( opt, &t8 );
        for( auto& c: json )
        {
            std::size_t const n = p.write_some( true, &c, 1, jec );
            BOOST_TEST( !jec.failed() );
            BOOST_TEST( n == 1 );
        }
        p.write_some(false, nullptr, 0, jec);
        BOOST_TEST( !jec.failed() );
        BOOST_TEST( t1 == t8 );
#if defined(__GNUC__) && __GNUC__ < 5
# pragma GCC diagnostic pop
#endif
    }

    template<class T>
    void testParseInto( T const& t )
    {
        testParseIntoValue<T>( value_from(t) );
    }

    template<class T>
    void
    testParseIntoErrors( error e, value const& sample )
    {
#if defined(__GNUC__) && __GNUC__ < 5
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
        system::error_code ec;
        T t1{};
        std::string json = serialize(sample);

        parse_into<T>(t1, json, ec);
        BOOST_TEST( ec.failed() );
        BOOST_TEST( ec.has_location() );
        BOOST_TEST( ec == e );

        T t2{};
        ec = {};
        parser_for<T> p( parse_options{}, &t2 );
        for( auto& c: json )
        {
            std::size_t const n = p.write_some( true, &c, 1, ec );
            if( ec.failed() )
                break;
            BOOST_TEST( n == 1 );
        }
        if( !ec.failed() )
            p.write_some(false, nullptr, 0, ec);

        BOOST_TEST( ec.failed() );
        BOOST_TEST( ec.has_location() );
        BOOST_TEST( ec == e );
#if defined(__GNUC__) && __GNUC__ < 5
# pragma GCC diagnostic pop
#endif
    }

    void testNull()
    {
        testParseInto( nullptr );
        testParseIntoErrors< std::nullptr_t >( error::not_null, 100 );
    }

    void testBoolean()
    {
        testParseInto( false );
        testParseInto( true );
        testParseIntoErrors< bool >( error::not_bool, nullptr );
    }

    void testIntegral()
    {
        testParseInto<char>( 'A' ); // ?
        testParseInto<signed char>( -127 );
        testParseInto<unsigned char>( 255 );
        testParseInto<short>( -32767 );
        testParseInto<unsigned short>( 65535 );
        testParseInto<int>( -32767 );
        testParseInto<unsigned int>( 65535 );
        testParseInto<long>( LONG_MIN );
        testParseInto<unsigned long>( ULONG_MAX );
        testParseInto<long long>( LLONG_MIN );
        testParseInto<unsigned long long>( ULLONG_MAX );

        testParseIntoErrors< int >( error::not_integer, "123" );
        testParseIntoErrors< int >( error::not_integer, true );
        testParseIntoErrors< int >( error::not_exact, LLONG_MIN );
        testParseIntoErrors< int >( error::not_exact, ULONG_MAX );
    }

    void testFloatingPoint()
    {
        testParseInto( 0.25f );
        testParseInto( 1.125 );
        // value_from doesn't support long double
        // testParseInto( 2.25L );
        {
            double d1 = 12;
            std::string json = serialize( value_from( 12 ) );

            system::error_code ec;
            double d2;
            parse_into(d2, json, ec);
            BOOST_TEST( !ec.failed() );
            BOOST_TEST( d1 == d2 );

            d1 = double(UINT64_MAX);
            json = serialize( value_from( UINT64_MAX ) );
            parse_into(d2, json, ec);
            BOOST_TEST( !ec.failed() );
            BOOST_TEST( d1 == d2 );
        }

        testParseIntoErrors< double >( error::not_double, { {"value", 12.1 } } );
    }

    void testString()
    {
        testParseInto<std::string>( "" );
        testParseInto<std::string>( "12345" );

        testParseIntoErrors< string >( error::not_string, UINT64_MAX );

        std::string s;
        parse_into(s, R"( "foobar" )");
        BOOST_TEST( s.size() == 6 );

        parse_into(s, R"( "qwertyuiop" )");
        BOOST_TEST( s.size() == 10 );
    }

    void testSequence()
    {
        testParseInto<std::vector<std::nullptr_t>>( { nullptr, nullptr } );

        testParseInto< std::vector<bool> >( {} );
        testParseInto< std::vector<bool> >( { true, false } );

        testParseInto< std::vector<int> >( {} );
        testParseInto< std::vector<int> >( { 1, 2, 3 } );
        testParseInto< std::vector<std::uint64_t> >( { 1, 2, UINT64_MAX } );

        testParseInto< std::vector<float> >( {} );
        testParseInto< std::vector<float> >( { 1.02f, 2.11f, 3.14f } );

        testParseInto< std::vector<std::string> >( {} );
        testParseInto< std::vector<std::string> >( { "one", "two", "three" } );

        testParseInto< std::vector<std::vector<int>> >( {} );
        testParseInto< std::vector<std::vector<int>> >( { {}, { 1 }, { 2, 3 }, { 4, 5, 6 } } );

        testParseInto< std::vector<std::map<std::string, int>> >(
            { {}, { {"1", 2}, {"3", 4} }, { { "5", 6 } } } );

        // clang <= 5 doesn't like when std::array is created from init-list
        std::array<int, 4> arr;
        arr.fill(17);
        testParseInto< std::array<int, 4> >( arr );

        testParseIntoErrors< std::vector<int> >( error::not_array, 1 );
        testParseIntoErrors< std::vector<char> >( error::not_array, "abcd" );
        testParseIntoErrors< std::vector<std::vector<int>> >(
            error::not_array, {1, 2, 3} );
        testParseIntoErrors< std::array<int, 4> >(
            error::size_mismatch, {1, 2, 3} );
        testParseIntoErrors< std::array<int, 4> >(
            error::size_mismatch, {1, 2, 3, 4, 5, 6, 7, 8} );

        testParseInto< std::vector<std::array<int, 4>> >( {arr,arr,arr} );

        std::vector<int> v;
        parse_into(v, "[1,2,3,4]");
        BOOST_TEST( v.size() == 4 );

        parse_into(v, "[5,6,7]");
        BOOST_TEST( v.size() == 3 );
    }

    void testMap()
    {
        testParseInto< std::map<std::string, int> >( {} );
        testParseInto< std::map<std::string, int> >( { { "one", 1 }, { "two", 2 } } );

        testParseInto< std::map<std::string, int> >( {} );
        testParseInto< std::map<std::string, int> >( { { "one", 1 }, { "two", 2 } } );

        testParseInto< std::map<std::string, std::vector<int>> >( {} );
        testParseInto< std::map<std::string, std::vector<int>> >( { { "one", { 1 } }, { "two", { 2, 3 } } } );

        testParseInto< std::map<std::string, std::map<std::string, int>> >( {} );
        testParseInto< std::map<std::string, std::map<std::string, int>> >( { { "one", {} }, { "two", { { "1", 1 }, { "2", 2 } } } } );

        testParseIntoErrors< std::map<std::string, int> >(
            error::not_object, { "1", 1, "2", 2} );
        testParseIntoErrors< std::map<std::string, std::map<std::string, int>> >(
            error::not_object, { {"1", {}}, {"2", {"3", 4}} } );

        std::map<std::string, int> m;
        parse_into(m, R"( {"1": 1, "2": 2, "3": 3} )");
        BOOST_TEST( m.size() == 3 );

        parse_into(m, R"( {"4": 4, "5": 5} )");
        BOOST_TEST( m.size() == 2 );
    }

    void testTuple()
    {
        testParseInto<std::pair<std::nullptr_t, std::uint64_t>>(
            {nullptr, UINT64_MAX} );

        testParseInto<std::tuple<int, float, std::string>>( {} );
        testParseInto<std::tuple<int, float, std::string>>( std::make_tuple(1, 3.14f, "hello") );

        testParseInto<std::vector<std::pair<int, int>>>( {} );
        testParseInto<std::vector<std::pair<int, int>>>( { { 1, 2 }, { 3, 4 } } );

        testParseInto<std::vector<std::vector<std::pair<int, int>>>>( {} );
        testParseInto<std::vector<std::vector<std::pair<int, int>>>>( { { { 1, 2 }, { 3, 4 } } } );
        testParseInto<std::vector<std::vector<std::pair<int, int>>>>( { { { 1, 2 }, { 3, 4 } }, { { 5, 6 }, { 7, 8 } } } );

        testParseInto<std::map<std::string, std::vector<std::pair<int, int>>>>( {} );
        testParseInto<std::map<std::string, std::vector<std::pair<int, int>>>>( { { "one", {} } } );
        testParseInto<std::map<std::string, std::vector<std::pair<int, int>>>>( { { "one", { { 1, 2 }, { 3, 4 } } } } );

        testParseInto<std::pair<std::vector<int>, std::map<std::string, std::pair<int, bool>>>>( {} );
        testParseInto<std::pair<std::vector<int>, std::map<std::string, std::pair<int, bool>>>>( { { 1, 2, 3 }, { { "one", { 7, true } } } } );

        testParseIntoErrors< std::pair<int, int> >( error::not_array, 1 );
        testParseIntoErrors< std::pair<int, int> >(
            error::size_mismatch, {1, 2, 3} );
        testParseIntoErrors< std::tuple<int, int, int> >(
            error::size_mismatch, {1, 2} );
        testParseIntoErrors< std::tuple<std::vector<int>>  >(
            error::size_mismatch, {{1,2}, {3,4}} );
        testParseIntoErrors<std::map<std::string, std::tuple<int, int>>>(
            error::size_mismatch, { {"tup", array()} });
    }

    void testStruct()
    {
#if defined(BOOST_DESCRIBE_CXX14)
        testParseInto<X>( {} );
        testParseInto<X>( { 1, 3.14f, "hello" } );

        testParseInto<Y>( {} );
        testParseInto<Y>( { { { 1, 1.0f, "one" }, { 2, 2.0f, "two" } }, { { "one", { 1, 1.1f, "1" } }, { "two", { 2, 2.2f, "2" } } } } );

        testParseInto<Z>( {} );
        testParseInto<Z>( { {1, 3.14f, "hello"}, true } );

        testParseIntoErrors<X>( error::not_object, 1 );
        testParseIntoErrors<X>( error::size_mismatch, { {"a", 1} } );

        object jo{ {"a", 1}, {"b", 3.14f}, {"c", "hello"} };
        jo["e1"] = array{1, ULONG_MAX, "three", array{}, true, .5, nullptr};
        jo["e2"] = object{ {"one", 1}, {"two", 2}, {"three", object{}} };
        jo["e3"] = "third extra element";
        jo["e4"] = 100;
        jo["e5"] = ULONG_MAX;
        jo["e6"] = 12.4;
        jo["e7"] = false;
        jo["e8"] = nullptr;
        testParseIntoValue<X>(jo);
#endif
    }

    void testEnum()
    {
#ifdef BOOST_DESCRIBE_CXX14
        testParseInto<E>( E::x );
        testParseInto<E>( E::y );
        testParseInto<E>( E::z );

        testParseIntoErrors< E >( error::not_string, (int)(E::y) );
        testParseIntoErrors< E >( error::unknown_name, "zoom" );
#endif // BOOST_DESCRIBE_CXX14
    }

    template< template<class...> class Variant, class Monostate >
    void testVariant()
    {
        testParseInto< Variant<int> >( 1 );
        testParseInto< Variant<int, std::string> >( 1 );
        testParseInto< Variant<int, std::string> >( "qwerty" );
        testParseInto< Variant<Monostate, int, std::string> >( {} );
        testParseInto< Variant<Monostate, int, std::string> >( 1 );
        testParseInto< Variant<Monostate, int, std::string> >( "qwerty" );
        testParseInto< Variant<bool, std::uint64_t> >( true );
        testParseInto< Variant<bool, std::uint64_t> >( UINT64_MAX );

        testParseInto< Variant< std::vector<int> > >(
            std::vector<int>{1, 2, 3, 4, 5} );
        testParseInto< Variant< Monostate, std::vector<int> > >(
            std::vector<int>{1, 2, 3, 4, 5} );

        testParseInto< std::vector< Variant<int, std::string> > >(
            {1, 2, 3, "four", 5, "six", "seven", 8});

        using V = Variant<
            std::vector< int >,
            std::tuple< int, std::string, std::map<std::string, int> >,
            std::tuple< int, std::string, std::map<std::string, double> > >;
        testParseInto< V >(
            std::make_tuple(
                5,
                "five",
                std::map<std::string, double>{ {"one", 1}, {"pi", 3.14} }));

        testParseIntoErrors< Variant<Monostate> >(
            error::exhausted_variants, "a" );
        testParseIntoErrors< Variant<int> >(
            error::exhausted_variants, "a" );
        testParseIntoErrors< Variant<Monostate, int, bool> >(
            error::exhausted_variants, "a" );
    }

    void testOptional()
    {
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
        testParseInto< std::optional<int> >( std::nullopt );
        testParseInto< std::optional<int> >( 1 );
        testParseInto< std::optional<std::uint64_t> >( ULONG_MAX );
        testParseInto< std::optional<bool> >( true );
        testParseInto< std::optional<double> >( 33.77 );

        testParseInto< std::optional<std::vector<std::nullptr_t>> >(
            std::nullopt );
        testParseInto< std::optional<std::vector<std::nullptr_t>> >(
           std::vector<std::nullptr_t>{} );
        testParseInto< std::optional<std::vector<std::nullptr_t>> >(
           std::vector<std::nullptr_t>{nullptr, nullptr} );

        testParseInto< std::optional<std::vector<std::string>> >(
           std::vector<std::string>{"1", "2", "3"} );

        testParseInto< std::optional< std::map<std::string, int> > >(
           std::map<std::string, int>{ {"1", 1}, {"2", 2}, {"3", 3} } );

        testParseInto< std::vector< std::optional<int> > >(
            {1, 2, 3, std::nullopt, 5, std::nullopt, std::nullopt, 8});
#endif
    }

    void testPath()
    {
#ifndef BOOST_NO_CXX17_HDR_FILESYSTEM
        testParseInto< std::filesystem::path >( "c:/" );
        testParseInto< std::filesystem::path >( ".." );
        testParseInto< std::filesystem::path >( "../" );
#endif // BOOST_NO_CXX17_HDR_FILESYSTEM
    }

    void run()
    {
        testNull();
        testBoolean();
        testIntegral();
        testFloatingPoint();
        testString();
        testSequence();
        testMap();
        testTuple();
        testStruct();
        testEnum();
        testOptional();
        testPath();
        testVariant<variant2::variant, variant2::monostate>();
#ifndef BOOST_NO_CXX17_HDR_VARIANT
        testVariant<std::variant, std::monostate>();
#endif

        {
            int n;
            BOOST_TEST_THROWS_WITH_LOCATION( parse_into( n, "null" ) );

            std::stringstream is("null");
            BOOST_TEST_THROWS_WITH_LOCATION( parse_into( n, is ) );
        }

        {
            int n;
            system::error_code ec;
            parse_into( n, "12 1", ec);
            BOOST_TEST( ec == error::extra_data );
            BOOST_TEST( ec.has_location() );
        }

        {
            std::stringstream is("12 1");
            int n;
            system::error_code ec;
            parse_into( n, is, ec);
            BOOST_TEST( ec == error::extra_data );
            BOOST_TEST( ec.has_location() );
        }

        {
            int n;
            std::stringstream is("1");
            is.setstate( std::ios::failbit );
            system::error_code ec;
            parse_into(n, is, ec);
            BOOST_TEST( ec == error::input_error );
            BOOST_TEST( ec.has_location() );
        }
    }
};

TEST_SUITE(parse_into_test, "boost.json.parse_into");

} // namespace boost
} // namespace json
