//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/value_to.hpp>

#include <boost/json/value_from.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/describe/class.hpp>
#include <boost/describe/enum.hpp>
#include <boost/variant2/variant.hpp>
#include <boost/config.hpp>

#include "test_suite.hpp"

#include <array>
#include <map>
#include <unordered_map>
#include <vector>

#ifndef BOOST_NO_CXX17_HDR_VARIANT
# include <variant>
#endif

#ifndef BOOST_NO_CXX17_HDR_FILESYSTEM
# include <filesystem>
#endif // BOOST_NO_CXX17_HDR_FILESYSTEM

namespace value_to_test_ns
{

struct T1 { };

//----------------------------------------------------------
struct T2 { };

boost::system::result<T2>
tag_invoke(
    boost::json::try_value_to_tag<T2>,
    boost::json::value const& jv)
{
    boost::json::string const* str = jv.if_string();
    if( str && *str == "T2" )
        return T2{};
    return make_error_code(boost::json::error::syntax);
}

//----------------------------------------------------------
struct T3 { };

T3
tag_invoke(
    boost::json::value_to_tag<T3>,
    boost::json::value const& jv)
{
    boost::json::string const* str = jv.if_string();
    if( !str )
        throw boost::system::system_error(
            make_error_code(boost::json::error::not_string));
    if ( *str != "T3" )
        throw std::invalid_argument("");
    return T3{};
}

// map-like type with fixed size
struct T4 {
    using value_type = std::pair<std::string, int>;

    value_type*
    begin()
    {
        return data;
    }

    value_type*
    end()
    {
        return data + sizeof(data);
    }

    std::pair< value_type*, bool >
    emplace(value_type);

    value_type data[2];
};

struct T5 { };

T5
tag_invoke(
    boost::json::value_to_tag<T5>,
    boost::json::value const& jv)
{
    // this is to shut up MSVC claiming
    // that this leads to unreachable code (duh)
    if( jv.is_object() && jv.get_object().size() > 1000000 )
        return T5{};

    throw std::bad_alloc();
}

//----------------------------------------------------------

struct T6
{
    int n;
    double d;
};
BOOST_DESCRIBE_STRUCT(T6, (), (n, d))

//----------------------------------------------------------

struct T7 : T6
{
    std::string s;

    bool
    get_b() const
    {
        return b;
    }

private:
    bool b = false;

    BOOST_DESCRIBE_CLASS(T7, (T6), (s), (), (b))
};

//----------------------------------------------------------

struct T10
{
    int n;
    T3 t3;
};
BOOST_DESCRIBE_STRUCT(T10, (), (n, t3))

//----------------------------------------------------------

BOOST_DEFINE_ENUM_CLASS(E1, a, b, c)

//----------------------------------------------------------

struct T8
{
    int n;
    double d;
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
    std::optional<std::string> opt_s;
#else
    std::string opt_s;
#endif // BOOST_NO_CXX17_HDR_OPTIONAL
};
BOOST_DESCRIBE_STRUCT(T8, (), (n, d, opt_s))

//----------------------------------------------------------

struct custom_context
{ };

struct T9
{ };

boost::system::result<T9>
tag_invoke(
    boost::json::try_value_to_tag<T9>,
    boost::json::value const& jv,
    custom_context const&)
{
    boost::json::string const* str = jv.if_string();
    if( str && *str == "T9" )
        return T9{};
    return make_error_code(boost::json::error::syntax);
}

// not default-constructible
struct T11
{
    int n;

    explicit T11( int v )
        : n(v)
    {}
};
T11
tag_invoke(
    boost::json::value_to_tag<T11>,
    boost::json::value const& jv)
{
    return T11( jv.to_number<int>() );
}

} // namespace value_to_test_ns

namespace std
{

// some versions of libstdc++ forward-declare tuple_size as class
#if defined(__clang__) || ( defined(__GNUC__) && __GNUC__ >= 10 )
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmismatched-tags"
#endif
template<>
struct tuple_size<value_to_test_ns::T4>
    : std::integral_constant<std::size_t, 2>
{ };
#if defined(__clang__) || ( defined(__GNUC__) && __GNUC__ >= 10 )
# pragma GCC diagnostic pop
#endif

} // namespace std

namespace boost {
namespace json {

template<>
struct is_null_like<::value_to_test_ns::T1> : std::true_type { };

template<>
struct is_described_class<::value_to_test_ns::T7> : std::true_type { };

template <class T, class = void>
struct can_apply_value_to
    : std::false_type
{
};

template <class T>
struct can_apply_value_to<T, detail::void_t<decltype(
    value_to<int>(std::declval<T>()))
>>
    : std::true_type
{
};

BOOST_STATIC_ASSERT(!can_apply_value_to<int>::value);

class value_to_test
{
public:

#define BOOST_TEST_CONV(x, ... ) \
    BOOST_TEST( value_to<decltype(x)>( value_from((x)), __VA_ARGS__ ) == (x) )

    template< class... Context >
    static
    void
    testNumberCast( Context const& ... ctx )
    {
        BOOST_TEST_CONV( (short)-1, ctx... );
        BOOST_TEST_CONV( (int)-2, ctx... );
        BOOST_TEST_CONV( (long)-3, ctx... );
        BOOST_TEST_CONV( (long long)-4, ctx... );
        BOOST_TEST_CONV( (unsigned short)1, ctx... );
        BOOST_TEST_CONV( (unsigned int)2, ctx... );
        BOOST_TEST_CONV( (unsigned long)3, ctx... );
        BOOST_TEST_CONV( (unsigned long long)4, ctx... );
        BOOST_TEST_CONV( (float)1.5, ctx... );
        BOOST_TEST_CONV( (double)2.5, ctx... );
        BOOST_TEST_CONV( true, ctx... );
        BOOST_TEST_THROWS_WITH_LOCATION(value_to<bool>(
            value(), ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(value_to<int>(
            value(), ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(value_to<double>(
            value(), ctx... ));
    }

    template< class... Context >
    static
    void
    testJsonTypes( Context const& ... ctx )
    {
        value_to<object>( value(object_kind), ctx... );
        value_to<array>( value(array_kind), ctx... );
        value_to<string>( value(string_kind), ctx... );

        BOOST_TEST_THROWS_WITH_LOCATION(value_to<object>(
            value(), ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(value_to<array>(
            value(), ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(value_to<string>(
            value(), ctx... ));
    }

    template< class... Context >
    static
    void
    testGenerics( Context const& ... ctx )
    {
        BOOST_TEST_CONV( std::string("test"), ctx... );
        BOOST_TEST_CONV(
            (std::map<std::string, int>
            {
                {"a", 1}, {"b", 2}, {"c", 3}
            }),
            ctx... );
        BOOST_TEST_CONV(
            (std::multimap<std::string, int>
            {
                {"2", 4}, {"3", 9}, {"5", 25}
            }),
            ctx... );
        BOOST_TEST_CONV(
            (std::unordered_map<std::string, int>
            {
                { "a", 1 }, {"b", 2}, {"c", 3}
            }),
            ctx... );
        BOOST_TEST_CONV( (std::vector<int>{1, 2, 3, 4}), ctx... );
        BOOST_TEST_CONV(
            (std::vector<bool>{true, false, false, true}), ctx... );
        BOOST_TEST_CONV(
            std::make_pair(std::string("test"), 5), ctx... );
        BOOST_TEST_CONV(
            std::make_tuple(
                std::string("outer"), std::make_pair(std::string("test"), 5) ),
            ctx... );
        BOOST_TEST_CONV(
            (std::map<int, int>
            {
                {2, 4}, {3, 9}, {5, 25}
            }),
            ctx... );

        {
            std::array<int, 500> arr;
            arr.fill(0);
            BOOST_TEST_CONV( arr, ctx... );
        }

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::string>( value(), ctx... ));

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<std::map<std::string, int>>( value(), ctx... )));
        // element fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<std::map<std::string, int>>(
                value{{"1", 1}, {"2", true}, {"3", false}}, ctx... )));
        // reserve fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<value_to_test_ns::T4>(
                value{{"1", 1}, {"2", true}, {"3", false}}, ctx... )));

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::vector<int>>( value(), ctx... ));
        // element fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::vector<int>>( value{1, 2, false, 3}, ctx... ));
        // reserve fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<std::array<int, 4>>( value{1, 2, 3}, ctx... )));
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<std::array<int, 4>>(
                value{1, 2, 3, 4, 5}, ctx... )));

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<std::tuple<int, int, int, int>>(
                value(), ctx... )));
        // element fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<std::tuple<int, int, int>>(
                value{1, 2, false}, ctx... )));
        // reserve fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            (value_to<std::tuple<int, int, int, int>>(
                value{1, 2, 3}, ctx... )));

    }

    void
    testContainerHelpers()
    {
        {
            std::vector<int> v;
            detail::try_reserve(
                v, 10, detail::reserve_implementation<decltype(v)>());
            BOOST_TEST(v.capacity() >= 10);
            BOOST_STATIC_ASSERT(std::is_same<
                decltype(detail::inserter(
                    v, detail::inserter_implementation<decltype(v)>())),
                decltype(std::back_inserter(v)) >::value);
        }
        {
            std::array<int, 2> arr;
            detail::try_reserve(
                arr, 2, detail::reserve_implementation<decltype(arr)>());
        }
        {
            int n;
            detail::try_reserve(
                n, 5, detail::reserve_implementation<decltype(n)>());
        }
    }

    template< class... Context >
    static
    void testNullptr( Context const& ... ctx )
    {
        (void)value_to<std::nullptr_t>( value(), ctx... );
        (void)value_to<::value_to_test_ns::T1>( value(), ctx... );
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::nullptr_t>( value(1), ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T1>( value(1), ctx... ));
    }

    template< class... Context >
    static
    void testDescribed( Context const& ... ctx )
    {
        ignore_unused( ctx... );
#ifdef BOOST_DESCRIBE_CXX14
        {
            value jv = {{"n", -78}, {"d", 0.125}};
            auto res = try_value_to<::value_to_test_ns::T6>(
                jv, ctx... );
            BOOST_TEST( res );
            BOOST_TEST( res->n == -78 );
            BOOST_TEST( res->d == 0.125 );

            jv.as_object()["x"] = 0;
            res = try_value_to<::value_to_test_ns::T6>(
                jv, ctx... );
            BOOST_TEST( res );
        }
        {
            value jv = {{"n", 1}, {"d", 2}, {"s", "xyz"}, {"b", true}};
            auto res = try_value_to<::value_to_test_ns::T7>(
                jv, ctx... );
            BOOST_TEST( res );
            BOOST_TEST( res->n == 1 );
            BOOST_TEST( res->d == 2 );
            BOOST_TEST( res->s == "xyz" );
            BOOST_TEST( res->get_b() == true );
        }

        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T6>( value(1), ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T6>( value{{"x", 0}}, ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T6>(
                value{{"n", 0}, {"x", 0}}, ctx... ));

        {
            value jv = "a";
            auto e1 = value_to<::value_to_test_ns::E1>( jv, ctx... );
            BOOST_TEST( e1 == ::value_to_test_ns::E1::a );

            jv = "b";
            e1 = value_to<::value_to_test_ns::E1>( jv, ctx... );
            BOOST_TEST( e1 == ::value_to_test_ns::E1::b );
        }

        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::E1>( value(1), ctx... ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::E1>( value("x"), ctx... ));

        {
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
            value jv = {{"n", -78}, {"d", 0.125}};
            auto res = try_value_to<::value_to_test_ns::T8>(
                jv, ctx... );
            BOOST_TEST( res );
            BOOST_TEST( res->n == -78 );
            BOOST_TEST( res->d == 0.125 );
            BOOST_TEST( std::nullopt == res->opt_s );

            jv.as_object()["x"] = 0;
            res = try_value_to<::value_to_test_ns::T8>(
                jv, ctx... );
            BOOST_TEST( res );
#endif // BOOST_NO_CXX17_HDR_OPTIONAL
        }

        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T10>(
                value{{"n", 0}, {"t3", "t10"}}, ctx... ));
#endif // BOOST_DESCRIBE_CXX14
    }

    template< class... Context >
    static
    void testOptional( Context const& ... ctx )
    {
        ignore_unused( ctx... );
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
        using Opts = std::vector<std::optional<int>>;
        value jv = value{1, nullptr, 3, nullptr, 5};
        auto opts = value_to<Opts>( jv, ctx... );
        BOOST_TEST( opts == (Opts{1, {}, 3, {}, 5}) );

        value_to< std::nullopt_t >(value());
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to< std::nullopt_t >( jv, ctx... ));
#endif
    }

    template<
        template <class...> class Variant, class Monostate, class... Context >
    static
    void
    testVariant( Context const& ... ctx )
    {
        using std::get;
        using Var = Variant<int, ::value_to_test_ns::T2, std::string>;

        value jv(4);
        auto v = value_to<Var>( jv, ctx... );
        BOOST_TEST( v.index() == 0 );
        BOOST_TEST( get<0>(v) == 4 );

        jv = "foobar";
        v = value_to<Var>( jv, ctx... );
        BOOST_TEST( v.index() == 2 );
        BOOST_TEST( get<2>(v) == "foobar" );

        jv = "T2";
        v = value_to<Var>( jv, ctx... );
        BOOST_TEST( v.index() == 1 );

        jv = 3.5;
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<Var>( jv, ctx... ));

        value_to<Monostate>( value(), ctx... );
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<Monostate>( jv, ctx... ));

        jv = 1024;
        using VT11 = Variant< value_to_test_ns::T11 >;
        VT11 v11 = value_to< VT11 >( jv, ctx... );
        BOOST_TEST( v11.index() == 0 );
        BOOST_TEST( get<0>(v11).n == 1024 );

        jv = nullptr;
        using V_T3_T1 = Variant<value_to_test_ns::T3, value_to_test_ns::T1>;
        auto v_t3_t1 = value_to<V_T3_T1>( jv, ctx... );
        BOOST_TEST( v_t3_t1.index() == 1 );
    }

    template< class... Context >
    static
    void testPath( Context const& ... ctx )
    {
        ignore_unused( ctx... );
#ifndef BOOST_NO_CXX17_HDR_FILESYSTEM
        using Paths = std::vector<std::filesystem::path>;
        value jv = value{"from/here", "to/there", "", "c:/" , "..", "../"};
        auto paths = value_to<Paths>( jv, ctx... );
        BOOST_TEST(
            paths == (Paths{
                "from/here", "to/there", "", "c:/" , "..", "../"}) );
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::filesystem::path>( value(1), ctx... ));
#endif // BOOST_NO_CXX17_HDR_FILESYSTEM
    }

    template< class... Context >
    static
    void
    testNonThrowing( Context const& ... ctx )
    {
        // using result
        {
            // clang 3.8 seems to have some bug when dealing with a lot of
            // template instantiations; this assert magically makes the problem
            // go away, I assume, by instantiating the needed types beforehand
            BOOST_STATIC_ASSERT(
                detail::conversion_round_trips<
                    mp11::mp_first<
                        mp11::mp_list<
                            Context..., int> >,
                    ::value_to_test_ns::T2,
                    detail::value_to_conversion>::value );

            auto res = try_value_to<::value_to_test_ns::T2>(
                value(), ctx... );
            BOOST_TEST( res.has_error() );
            BOOST_TEST( res.error() == error::syntax );

            res = try_value_to<::value_to_test_ns::T2>(
                value("T2"), ctx... );
            BOOST_TEST( res.has_value() );
        }
        // throwing overload falls back to nonthrowing customization
        {
            BOOST_TEST_THROWS(
                value_to<::value_to_test_ns::T2>( value(), ctx... ),
                system::system_error);
        }
        // nonthrowing overload falls back to throwing customization
        {
            auto res = try_value_to<::value_to_test_ns::T3>(
                value(), ctx... );
            BOOST_TEST( res.has_error() );
            BOOST_TEST( res.error() == error::not_string );

            res = try_value_to<::value_to_test_ns::T3>(
                value(""), ctx... );
            BOOST_TEST( res.has_error() );
            BOOST_TEST( res.error() == error::exception );

            res = try_value_to<::value_to_test_ns::T3>(
                value("T3"), ctx... );
            BOOST_TEST( res.has_value() );
        }
        // sequence
        {
            // wrong input type
            {
                auto res = try_value_to< std::vector<::value_to_test_ns::T2> >(
                    value("not an array"), ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::not_array );
            }

            // wrong input type
            {
                auto res = try_value_to< std::array<int, 4> >(
                    value{1, 2}, ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::size_mismatch );
            }

            // element error
            {
                auto res = try_value_to< std::vector<::value_to_test_ns::T2> >(
                    value{"T2", "T2", nullptr}, ctx... );
                BOOST_TEST( res.error() == error::syntax );
            }

            // success
            auto res = try_value_to< std::vector<::value_to_test_ns::T2> >(
                value{"T2", "T2", "T2"}, ctx... );
            BOOST_TEST( res.has_value() );
            BOOST_TEST( res->size() == 3 );
        }
        // map
        {
            // wrong input type
            {
                auto res = try_value_to<
                    std::map<std::string, ::value_to_test_ns::T2> >(
                        value("not a map"), ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::not_object );
            }

            // reserve fails
            {
                auto res = try_value_to< ::value_to_test_ns::T4 >(
                    value{{"1", 1}, {"2", 2}, {"3", 3}}, ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::size_mismatch );
            }

            // element error
            {
                auto res = try_value_to<
                    std::map<std::string, ::value_to_test_ns::T2> >(
                        value{{"1", "T2"}, {"2", "T2"}, {"3", nullptr}},
                        ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error() == error::syntax );
            }

            // success
            auto res = try_value_to<
                std::map<std::string, ::value_to_test_ns::T2> >(
                    value{{"1", "T2"}, {"2", "T2"}, {"3", "T2"}}, ctx... );
            BOOST_TEST( res.has_value() );
            BOOST_TEST( res->size() == 3 );
        }
        // tuple
        {
            // wrong input type
            {
                auto res = try_value_to<std::tuple<
                        int,
                        std::string,
                        bool,
                        std::nullptr_t,
                        ::value_to_test_ns::T2>>(
                    value("not an array"), ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::not_array );
            }

            // size mismatch
            {
                auto res = try_value_to<std::tuple<
                        int,
                        std::string,
                        bool,
                        std::nullptr_t,
                        ::value_to_test_ns::T2>>(
                    value{1, 2}, ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::size_mismatch );
            }

            // element error
            {
                auto res = try_value_to<std::tuple<
                        int,
                        std::string,
                        bool,
                        std::nullptr_t,
                        ::value_to_test_ns::T2>>(
                    value{1, "foobar", false, nullptr, ""}, ctx... );
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error() == error::syntax );
            }

            // success
            auto res = try_value_to<std::tuple<
                    int,
                    std::string,
                    bool,
                    std::nullptr_t,
                    ::value_to_test_ns::T2>>(
                value{1, "foobar", false, nullptr, "T2"}, ctx... );
            BOOST_TEST( res.has_value() );
            BOOST_TEST( std::get<0>(*res) == 1 );
            BOOST_TEST( std::get<1>(*res) == "foobar" );
            BOOST_TEST_NOT( std::get<2>(*res) );
        }
        // rethrowing bad_alloc
        BOOST_TEST_THROWS(
            try_value_to<value_to_test_ns::T5>( value(), ctx... ),
            std::bad_alloc);
    }

    template< class... Context >
    static
    void
    testUserConversion( Context const& ... ctx )
    {
        value_to<value_to_test_ns::T2>( value("T2"), ctx... );
    }

    void
    testContext()
    {
        value_to<value_to_test_ns::T9>(
            value("T9"), value_to_test_ns::custom_context() );

        BOOST_TEST_THROWS(
            value_to<value_to_test_ns::T9>(
                value(), value_to_test_ns::custom_context() ),
            system::system_error);
    }

    struct run_templated_tests
    {
        // this overload supports zero or one default constructible contexts
        // and used with mp_for_each
        template< class... Context >
        void operator()( mp11::mp_list< Context... > )
        {
            testNumberCast( Context()... );
            testJsonTypes( Context()... );
            testGenerics( Context()... );
            testNullptr( Context()... );
            testDescribed( Context()... );
            testOptional( Context()... );
            testVariant< variant2::variant, variant2::monostate > ( Context()... );
#ifndef BOOST_NO_CXX17_HDR_VARIANT
            testVariant< std::variant, std::monostate > ( Context()... );
#endif // BOOST_NO_CXX17_HDR_VARIANT
            testPath( Context()... );
            testNonThrowing( Context()... );
            testUserConversion( Context()... );
        }
    };

    void
    run()
    {
        mp11::mp_for_each<
            mp11::mp_list<
                mp11::mp_list<>,
                mp11::mp_list<detail::no_context>,
                mp11::mp_list<value_to_test_ns::custom_context>,
                mp11::mp_list<
                    std::tuple<value_to_test_ns::custom_context>>,
                mp11::mp_list<
                    std::tuple<
                        std::tuple<value_to_test_ns::custom_context>>>,
                mp11::mp_list<
                    std::tuple<
                        detail::no_context, value_to_test_ns::custom_context>>
             >>( run_templated_tests() );

        testContext();
        testContainerHelpers();
    }
};

TEST_SUITE(value_to_test, "boost.json.value_to");

} // namespace json
} // namespace boost
