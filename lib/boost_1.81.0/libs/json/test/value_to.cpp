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
#include <boost/describe/class.hpp>
#include <boost/describe/enum.hpp>

#include "test_suite.hpp"

#include <array>
#include <map>
#include <unordered_map>
#include <vector>

namespace value_to_test_ns
{

struct T1 { };

//----------------------------------------------------------
struct T2 { };

boost::json::result<T2>
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
        throw boost::json::system_error(
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
    boost::json::value const&)
{
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
};
BOOST_DESCRIBE_STRUCT(T7, (T6), (s))

//----------------------------------------------------------

BOOST_DEFINE_ENUM_CLASS(E1, a, b, c)

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

BOOST_JSON_NS_BEGIN

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

    template<class T>
    void
    check(T t)
    {
        BOOST_TEST(value_to<T>(value_from(t)) == t);
    }

    void
    testNumberCast()
    {
        check((short)-1);
        check((int)-2);
        check((long)-3);
        check((long long)-4);
        check((unsigned short)1);
        check((unsigned int)2);
        check((unsigned long)3);
        check((unsigned long long)4);
        check((float)1.5);
        check((double)2.5);
        check(true);
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<bool>(value()) );
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<int>(value()) );
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<double>(value()) );
    }

    void
    testJsonTypes()
    {
        value_to<object>(value(object_kind));
        value_to<array>(value(array_kind));
        value_to<string>(value(string_kind));

        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<object>(value()) );
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<array>(value()) );
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<string>(value()) );
    }

    void
    testGenerics()
    {
        check(std::string("test"));
        check(std::map<std::string, int>
        {
            {"a", 1}, {"b", 2}, {"c", 3}
        });
        check(std::multimap<std::string, int>
        {
            {"2", 4}, {"3", 9}, {"5", 25}
        });
        check(std::unordered_map<std::string, int>
        {
            { "a", 1 }, {"b", 2}, {"c", 3}
        });
        check(std::vector<int>{1, 2, 3, 4});
        check(std::vector<bool>{true, false, false, true});
        check(std::make_pair(std::string("test"), 5));
        check(std::make_tuple(std::string("outer"),
            std::make_pair(std::string("test"), 5)));
        check(std::map<int, int>
        {
            {2, 4}, {3, 9}, {5, 25}
        });

        {
            std::array<int, 1000> arr;
            arr.fill(0);
            check(arr);
        }

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::string>(value()) );

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::map<std::string, int>>(value()) ));
        // element fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::map<std::string, int>>(
                value{{"1", 1}, {"2", true}, {"3", false}}) ));
        // reserve fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<value_to_test_ns::T4>(
                value{{"1", 1}, {"2", true}, {"3", false}}) ));

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::vector<int>>(value()) );
        // element fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<std::vector<int>>(value{1, 2, false, 3}) );
        // reserve fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::array<int, 4>>(value{1, 2, 3}) ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::array<int, 4>>(value{1, 2, 3, 4, 5}) ));

        // mismatched type
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::tuple<int, int, int, int>>(value()) ));
        // element fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::tuple<int, int, int>>(value{1, 2, false}) ));
        // reserve fails
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::tuple<int, int, int, int>>(value{1, 2, 3}) ));

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

    void testNullptr()
    {
        (void)value_to<std::nullptr_t>(value());
        (void)value_to<::value_to_test_ns::T1>(value());
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<std::nullptr_t>(value(1)) ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            ( value_to<::value_to_test_ns::T1>(value(1)) ));
    }

    void testDescribed()
    {
#ifdef BOOST_DESCRIBE_CXX14
        {
            value jv = {{"n", -78}, {"d", 0.125}};
            auto res = try_value_to<::value_to_test_ns::T6>(jv);
            BOOST_TEST( res );
            BOOST_TEST( res->n == -78 );
            BOOST_TEST( res->d == 0.125 );
        }
        {
            value jv = {{"n", 1}, {"d", 2}, {"s", "xyz"}};
            auto res = try_value_to<::value_to_test_ns::T7>(jv);
            BOOST_TEST( res );
            BOOST_TEST( res->n == 1 );
            BOOST_TEST( res->d == 2 );
            BOOST_TEST( res->s == "xyz" );
        }

        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T6>( value(1) ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T6>( value{{"x", 0}} ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::T6>( (value{{"n", 0}, {"x", 0}}) ));

        {
            value jv = "a";
            auto e1 = value_to<::value_to_test_ns::E1>(jv);
            BOOST_TEST( e1 == ::value_to_test_ns::E1::a );

            jv = "b";
            e1 = value_to<::value_to_test_ns::E1>(jv);
            BOOST_TEST( e1 == ::value_to_test_ns::E1::b );
        }

        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::E1>( value(1) ));
        BOOST_TEST_THROWS_WITH_LOCATION(
            value_to<::value_to_test_ns::E1>( value("x") ));
#endif
    }

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
    void testOptional()
    {
        using Opts = std::vector<std::optional<int>>;
        value jv = value{1, nullptr, 3, nullptr, 5};
        auto opts = value_to<Opts>(jv);
        BOOST_TEST( opts == (Opts{1, {}, 3, {}, 5}) );

        value_to< std::nullopt_t >(value());
        BOOST_TEST_THROWS_WITH_LOCATION( value_to< std::nullopt_t >(jv) );
    }
#endif

#ifndef BOOST_NO_CXX17_HDR_VARIANT
    void
    testVariant()
    {
        using Var = std::variant<int, ::value_to_test_ns::T2, std::string>;

        value jv(4);
        auto v = value_to<Var>(jv);
        BOOST_TEST( v.index() == 0 );
        BOOST_TEST( std::get<0>(v) == 4 );

        jv = "foobar";
        v = value_to<Var>(jv);
        BOOST_TEST( v.index() == 2 );
        BOOST_TEST( std::get<2>(v) == "foobar" );

        jv = "T2";
        v = value_to<Var>(jv);
        BOOST_TEST( v.index() == 1 );

        jv = 3.5;
        BOOST_TEST_THROWS_WITH_LOCATION( value_to<Var>(jv) );

        value_to<std::monostate>( value() );
        BOOST_TEST_THROWS_WITH_LOCATION( value_to<std::monostate>(jv) );
    }
#endif // BOOST_NO_CXX17_HDR_VARIANT

    void
    testNonThrowing()
    {
        // using result
        {
            auto res = try_value_to<::value_to_test_ns::T2>(value());
            BOOST_TEST( res.has_error() );
            BOOST_TEST( res.error() == error::syntax );

            res = try_value_to<::value_to_test_ns::T2>(value("T2"));
            BOOST_TEST( res.has_value() );
        }
        // throwing overload falls back to nonthrowing customization
        {
            BOOST_TEST_THROWS(
                value_to<::value_to_test_ns::T2>(value()),
                system_error);
        }
        // nonthrowing overload falls back to throwing customization
        {
            auto res = try_value_to<::value_to_test_ns::T3>(value());
            BOOST_TEST( res.has_error() );
            BOOST_TEST( res.error() == error::not_string );

            res = try_value_to<::value_to_test_ns::T3>(value(""));
            BOOST_TEST( res.has_error() );
            BOOST_TEST( res.error() == error::exception );

            res = try_value_to<::value_to_test_ns::T3>(value("T3"));
            BOOST_TEST( res.has_value() );
        }
        // sequence
        {
            // wrong input type
            {
                auto res = try_value_to< std::vector<::value_to_test_ns::T2> >(
                    value("not an array"));
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::not_array );
            }

            // wrong input type
            {
                auto res = try_value_to< std::array<int, 4> >(
                    value{1, 2});
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::size_mismatch );
            }

            // element error
            {
                auto res = try_value_to< std::vector<::value_to_test_ns::T2> >(
                    value{"T2", "T2", nullptr});
                BOOST_TEST( res.error() == error::syntax );
            }

            // success
            auto res = try_value_to< std::vector<::value_to_test_ns::T2> >(
                value{"T2", "T2", "T2"});
            BOOST_TEST( res.has_value() );
            BOOST_TEST( res->size() == 3 );
        }
        // map
        {
            // wrong input type
            {
                auto res = try_value_to<
                        std::map<std::string, ::value_to_test_ns::T2> >(
                    value("not a map"));
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::not_object );
            }

            // reserve fails
            {
                auto res = try_value_to< ::value_to_test_ns::T4 >(
                    value{{"1", 1}, {"2", 2}, {"3", 3}});
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error().has_location() );
                BOOST_TEST( res.error() == error::size_mismatch );
            }

            // element error
            {
                auto res = try_value_to<
                        std::map<std::string, ::value_to_test_ns::T2> >(
                    value{{"1", "T2"}, {"2", "T2"}, {"3", nullptr}});
                BOOST_TEST( res.has_error() );
                BOOST_TEST( res.error() == error::syntax );
            }

            // success
            auto res = try_value_to<
                    std::map<std::string, ::value_to_test_ns::T2> >(
                value{{"1", "T2"}, {"2", "T2"}, {"3", "T2"}});
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
                    value("not an array"));
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
                    value{1, 2});
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
                    value{1, "foobar", false, nullptr, ""});
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
                value{1, "foobar", false, nullptr, "T2"});
            BOOST_TEST( res.has_value() );
            BOOST_TEST( std::get<0>(*res) == 1 );
            BOOST_TEST( std::get<1>(*res) == "foobar" );
            BOOST_TEST_NOT( std::get<2>(*res) );
        }
        // rethrowing bad_alloc
        BOOST_TEST_THROWS(
            try_value_to<value_to_test_ns::T5>(value()),
            std::bad_alloc);
    }

    void
    testUserConversion()
    {
        value_to<value_to_test_ns::T2>(value("T2"));
    }

    void
    run()
    {
        testNumberCast();
        testJsonTypes();
        testGenerics();
        testContainerHelpers();
        testNullptr();
        testDescribed();
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
        testOptional();
#endif
#ifndef BOOST_NO_CXX17_HDR_VARIANT
        testVariant();
#endif // BOOST_NO_CXX17_HDR_VARIANT
        testUserConversion();
        testNonThrowing();
    }
};

TEST_SUITE(value_to_test, "boost.json.value_to");

BOOST_JSON_NS_END
