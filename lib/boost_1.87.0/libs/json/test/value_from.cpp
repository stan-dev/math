//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
// Copyright (c) 2020 Krystian Stasiowski (sdkrystian@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/value_from.hpp>

#include <boost/json/value.hpp> // prevent intellisense bugs
#include <boost/json/serialize.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/describe/class.hpp>
#include <boost/describe/enum.hpp>

#include "test_suite.hpp"

#include <array>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <unordered_map>

#ifndef BOOST_NO_CXX17_HDR_FILESYSTEM
# include <filesystem>
#endif // BOOST_NO_CXX17_HDR_FILESYSTEM

//----------------------------------------------------------

namespace value_from_test_ns
{

//----------------------------------------------------------

struct T1
{
    int i = 42;
};

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T1 const& t)
{
    jv = t.i;
}

//----------------------------------------------------------

// member function
// uses generic algorithms
struct T2
{
    std::vector<int> v;
    std::string s;

    T2()
        : v({1,2,3, 4})
        , s("test")
    {
    }
};

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T2 const& t)
{
    jv = { boost::json::value_from(t.v), boost::json::value_from(t.s) };
}

struct T3
{
};

BOOST_STATIC_ASSERT(! ::boost::json::has_value_from<T3>::value);

struct T4
{
    operator ::boost::json::value()
    {
        return nullptr;
    }
};

BOOST_STATIC_ASSERT(! ::boost::json::has_value_from<T4>::value);

//----------------------------------------------------------

struct T5 : std::vector<int>
{
    using std::vector<int>::vector;
};

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T5 const&)
{
    jv = "T5";
}


//----------------------------------------------------------

struct T6
{
};

inline
std::size_t
size(T6 const&)
{
    return 3;
}

struct T7 { };

//----------------------------------------------------------

struct T8
{
    bool error;
};

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    ::boost::system::error_code& ec,
    T8 const& t8)
{
    if( t8.error )
    {
        ec = ::boost::json::error::syntax;
        return;
    }

    jv = "T8";
}

//----------------------------------------------------------

struct T9
{
    int num;
};

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T9 const& t9)
{
    if( t9.num == 0 )
        throw std::invalid_argument("");
    if( t9.num < 0 )
        throw ::boost::system::system_error(
            make_error_code(::boost::json::error::syntax));

    jv = "T9";
}

//----------------------------------------------------------

struct T10
{
    int n;
    double d;
};
BOOST_DESCRIBE_STRUCT(T10, (), (n, d))

//----------------------------------------------------------

struct T11 : T10
{
    std::string s;

    void
    set_b(bool v)
    {
        b = v;
    }

private:
    bool b;

    BOOST_DESCRIBE_CLASS(T11, (T10), (s), (), (b))
};

//----------------------------------------------------------

BOOST_DEFINE_ENUM_CLASS(E1, a, b, c)

//----------------------------------------------------------

struct custom_context
{ };

struct T12
{ };

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T12 const&,
    custom_context const& )
{
    jv = "T12";
}

//----------------------------------------------------------

struct another_context
{ };

struct T13
{ };

void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    T13 const&,
    another_context const& )
{
    jv = "T13";
}

//----------------------------------------------------------

template< class T, class Ctx >
void
tag_invoke(
    ::boost::json::value_from_tag,
    ::boost::json::value& jv,
    std::vector< T > const& vec,
    another_context const&,
    Ctx const& ctx )
{
    boost::json::object& jo = jv.emplace_object();
    unsigned i = 0;
    for( auto&& item: vec )
    {
        jo.emplace(
            std::to_string(++i),
            boost::json::value_from( item, ctx, jo.storage() ) );
    }
}


} // namespace value_from_test_ns

template<class T>
static
void
check(
    ::boost::json::string_view s,
    T const& t)
{
    {
        auto const jv = value_from(t,
            ::boost::json::storage_ptr{});
        auto const js =
            ::boost::json::serialize(jv);
        BOOST_TEST(js == s);
    }
    {
        auto const jv =
            ::boost::json::value_from(t);
        auto const js =
            ::boost::json::serialize(jv);
        BOOST_TEST(js == s);
    }
}

namespace boost {
namespace json {

template<>
struct is_null_like<::value_from_test_ns::T7>
    : std::true_type
{ };

template<>
struct is_described_class<::value_from_test_ns::T11>
    : std::true_type
{ };

namespace {

template< class T, class... Context >
static
void
testValueCtor( T const& t, Context const& ... ctx )
{
    BOOST_TEST(
        serialize(value_from( t, ctx... )) == serialize(value(t)) );
}

} // namespace

// integral
BOOST_STATIC_ASSERT(has_value_from<int>::value);
BOOST_STATIC_ASSERT(has_value_from<int&>::value);
BOOST_STATIC_ASSERT(has_value_from<int&&>::value);
// array
BOOST_STATIC_ASSERT(has_value_from<int[4]>::value);
BOOST_STATIC_ASSERT(has_value_from<int(&)[4]>::value);
BOOST_STATIC_ASSERT(has_value_from<int(&&)[4]>::value);
// forward range
BOOST_STATIC_ASSERT(has_value_from<std::vector<int>>::value);
BOOST_STATIC_ASSERT(has_value_from<std::vector<int>&>::value);
BOOST_STATIC_ASSERT(has_value_from<std::vector<int>&&>::value);
// tuple-like
BOOST_STATIC_ASSERT(has_value_from<std::tuple<int, int>>::value);
BOOST_STATIC_ASSERT(has_value_from<std::tuple<int, int>&>::value);
BOOST_STATIC_ASSERT(has_value_from<std::tuple<int, int>&&>::value);
BOOST_STATIC_ASSERT(has_value_from<key_value_pair>::value);
BOOST_STATIC_ASSERT(has_value_from<key_value_pair&>::value);
BOOST_STATIC_ASSERT(has_value_from<key_value_pair&&>::value);

// object-like
BOOST_STATIC_ASSERT(has_value_from<std::map<string_view, int>>::value);

class value_from_test
{
public:
    template< class... Context >
    static
    void
    testValueCtors( Context const& ... ctx )
    {
        // value_from supports every value constructor

        testValueCtor<value, Context...>( value(), ctx... );

        char const* s = "5";
        testValueCtor<char const*, Context...>( s, ctx... );
    }

    template< class... Context >
    static
    void
    testGeneral( Context const& ... ctx )
    {
        {
            int a[4] = {1, 2, 3, 4};
            value b{1, 2, 3, 4};
            value c = value_from( a, ctx... );
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            std::tuple<int, string, int, bool> a{1, "2", 42, true};
            value b{1, "2", 42, true};
            value c = value_from( a, ctx... );
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            std::array<int, 500> a;
            a.fill(0);

            value b;
            array& b_arr = b.emplace_array();
            b_arr.insert(b_arr.end(), a.begin(), a.end());

            BOOST_TEST(value_from( a, ctx... ) == b);
        }
        {
            std::pair<int, string> a{1, string("2")};
            value b{1, "2"};
            value c = value_from( a, ctx... );
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            // ensures that this isn't parsed as a key value pair
            std::pair<string_view, int> a{"2", 1};
            value b{"2", 1};
            value c = value_from( a, ctx... );
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            key_value_pair a{"2", 1};
            value b{"2", 1};
            value c = value_from( a, ctx... );
            BOOST_TEST(c.is_array());
            BOOST_TEST(serialize(c) == serialize(b));
        }
        {
            ::value_from_test_ns::T7 a;
            value b = value_from( a, ctx... );
            BOOST_TEST(b.is_null());
        }
    }

    template< class... Context >
    static
    void
    testAssociative( Context const& ... ctx )
    {
        {
            std::map<string_view, int> a =
                {{"a", 1}, {"b", 2}, {"c", 3}};
            value b = {{"a", 1}, {"b", 2}, {"c", 3}};
            value c = value_from( a, ctx...  );
            BOOST_TEST(c.is_object());
            BOOST_TEST(a.size() == c.as_object().size());
            BOOST_TEST(b.as_object().size() == c.as_object().size());
        }
        {
            std::unordered_map<std::string, int> a =
               {{"a", 1}, {"b", 2}, {"c", 3}};
            value b = {{"a", 1}, {"b", 2}, {"c", 3}};
            value c = value_from( a, ctx... );
            BOOST_TEST(c.is_object());
            BOOST_TEST(a.size() == c.as_object().size());
            BOOST_TEST(b.as_object().size() == c.as_object().size());
        }
        {
            std::map<int, int> a =
                {{1, 1}, {2, 2}, {3, 3}};
            value b = {{1, 1}, {2, 2}, {3, 3}};
            value c = value_from( a, ctx... );
            BOOST_TEST(!c.is_object());
            BOOST_TEST(a.size() == c.as_array().size());
            BOOST_TEST(b.as_array().size() == c.as_array().size());
        }
        {
            std::unordered_map<int, int> a =
                {{1, 1}, {2, 2}, {3, 3}};
            value b = {{1, 1}, {2, 2}, {3, 3}};
            value c = value_from( a, ctx... );
            BOOST_TEST(!c.is_object());
            BOOST_TEST(a.size() == c.as_array().size());
            BOOST_TEST(b.as_array().size() == c.as_array().size());
        }
    }

    template< class... Context >
    static
    void
    testPreferUserCustomizations( Context const& ... ctx )
    {
        value_from_test_ns::T5 t5;
        BOOST_TEST((::boost::json::value_from( t5, ctx...  ) == "T5"));
    }

    void
    testTrySize()
    {
        {
            std::vector<int> v{1, 2, 3};
            using impl = detail::size_implementation<decltype(v)>;
            BOOST_TEST(detail::try_size(v, impl()) == 3);
        }
        {
            int arr[4] = {1, 2, 3, 4};
            using impl = detail::size_implementation<decltype(arr)>;
            BOOST_TEST(detail::try_size(arr, impl()) == 4);
        }
        {
            value_from_test_ns::T6 t;
            using impl = detail::size_implementation<decltype(t)>;
            BOOST_TEST(detail::try_size(t, impl()) == 3);
        }
        {
            value_from_test_ns::T1 t;
            using impl = detail::size_implementation<decltype(t)>;
            BOOST_TEST(detail::try_size(t, impl()) == 0);
        }
    }

    template< class... Context >
    static
    void
    testDescribed( Context const& ... ctx )
    {
        ignore_unused( ctx... );
#ifdef BOOST_DESCRIBE_CXX14
        ::value_from_test_ns::T10 t10{909, -1.45};
        value jv = value_from( t10, ctx... );
        BOOST_TEST(( jv == value{{"n", 909}, {"d", -1.45}} ));

        ::value_from_test_ns::T11 t11;
        t11.n = 67;
        t11.d = -.12;
        t11.s = "qwerty";
        t11.set_b(true);
        jv = value_from( t11, ctx... );
        BOOST_TEST(( jv == value{{"n", 67}, {"d", -.12}, {"s", "qwerty"}, {"b", true}} ));

        ::value_from_test_ns::E1 e1 = ::value_from_test_ns::E1::a;
        BOOST_TEST( value_from( e1, ctx... ) == "a" );

        e1 = ::value_from_test_ns::E1::b;
        BOOST_TEST( value_from( e1, ctx... ) == "b" );

        e1 = static_cast<::value_from_test_ns::E1>(1001);
        BOOST_TEST( value_from( e1, ctx... ) == 1001 );
#endif
    }

    template< class... Context >
    static
    void testOptional( Context const& ... ctx )
    {
        ignore_unused( ctx... );
#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
        std::vector<std::optional<int>> opts{1, 2, 3, {}, 5};
        value jv = value_from( opts, ctx... );
        BOOST_TEST( jv == (value{1, 2, 3, nullptr, 5}) );

        BOOST_TEST( value_from( std::nullopt, ctx... ).is_null() );
#endif
    }

    template< class... Context >
    static
    void
    testVariant( Context const& ... ctx )
    {
        ignore_unused( ctx... );
#ifndef BOOST_NO_CXX17_HDR_VARIANT
        std::variant<int, ::value_from_test_ns::T5, double> v = 4;
        value jv = value_from( v, ctx... );
        BOOST_TEST(jv == 4);

        v = 0.5;
        jv = value_from( v, ctx... );
        BOOST_TEST(jv == 0.5);

        v = ::value_from_test_ns::T5{};
        jv = value_from( v, ctx... );
        BOOST_TEST(jv == "T5");

        BOOST_TEST( value() == value_from( std::monostate(), ctx... ) );
#endif // BOOST_NO_CXX17_HDR_VARIANT
    }

    template< class... Context >
    static
    void testPath( Context const& ... ctx )
    {
        ignore_unused( ctx... );
#ifndef BOOST_NO_CXX17_HDR_FILESYSTEM
        std::vector<std::filesystem::path> paths{
            "from/here", "to/there", "", "c:/" , "..", "../"};
        value jv = value_from( paths, ctx... );
        BOOST_TEST(
            jv == (value{"from/here", "to/there", "", "c:/" , "..", "../"}) );
#endif // BOOST_NO_CXX17_HDR_FILESYSTEM
    }

    void
    testContext()
    {
        value jv = value_from(
            value_from_test_ns::T1{}, value_from_test_ns::custom_context() );
        BOOST_TEST( jv == 42 );

        jv = value_from(
            value_from_test_ns::T12{},
            value_from_test_ns::custom_context() );
        BOOST_TEST( jv == "T12" );

        jv = value_from(
            value_from_test_ns::T13{},
            std::make_tuple(
                value_from_test_ns::custom_context(),
                value_from_test_ns::another_context() ) );
        BOOST_TEST( jv == "T13" );

        jv = value_from(
            std::pair<value_from_test_ns::T12, value_from_test_ns::T13>(),
            std::make_tuple(
                value_from_test_ns::custom_context(),
                value_from_test_ns::another_context() ) );
        BOOST_TEST( jv == (array{ "T12", "T13" }) );


        auto ctx = std::make_tuple(
            value_from_test_ns::custom_context(),
            value_from_test_ns::another_context() );
        using Ctx = decltype(ctx);
        using Sup = detail::supported_context<
            Ctx,
            std::vector<value_from_test_ns::T12>,
            detail::value_from_conversion >;
        BOOST_STATIC_ASSERT( Sup::index::value == 1 );
        BOOST_STATIC_ASSERT(
            std::is_same<
                Sup::type,
                value_from_test_ns::another_context>::value );

        jv = value_from(
            std::vector<value_from_test_ns::T12>(2),
            std::make_tuple(
                value_from_test_ns::custom_context(),
                value_from_test_ns::another_context() ) );
        BOOST_TEST( jv == (object{ {"1", "T12"}, {"2", "T12"} }) );
    }

    struct run_templated_tests
    {
        // this overload supports zero or one default constructible contexts
        // and used with mp_for_each
        template< class... Context >
        void operator()( mp11::mp_list< Context... > )
        {
            testValueCtors( Context()... );
            testGeneral( Context()... );
            testAssociative( Context()... );
            testPreferUserCustomizations( Context()... );
            testDescribed( Context()... );
            testOptional( Context()... );
            testVariant( Context()... );
            testPath( Context()... );
        }
    };

    void
    run()
    {
        check("42", ::value_from_test_ns::T1{});
        check("[[1,2,3,4],\"test\"]", ::value_from_test_ns::T2{});

        mp11::mp_for_each<
            mp11::mp_list<
                mp11::mp_list<>,
                mp11::mp_list<detail::no_context>,
                mp11::mp_list<value_from_test_ns::custom_context>,
                mp11::mp_list<
                    std::tuple<value_from_test_ns::custom_context>>,
                mp11::mp_list<
                    std::tuple<
                        std::tuple<value_from_test_ns::custom_context>>>,
                mp11::mp_list<
                    std::tuple<
                        detail::no_context,
                        value_from_test_ns::custom_context>>
            >>( run_templated_tests() );

        testContext();
        testTrySize();
    }
};

TEST_SUITE(value_from_test, "boost.json.value_from");

} // namespace json
} // namespace boost
