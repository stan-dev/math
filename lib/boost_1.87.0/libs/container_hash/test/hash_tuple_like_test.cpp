// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__clang__)
# pragma clang diagnostic ignored "-Wmismatched-tags"
#endif

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/type_traits/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/config.hpp>
#include <utility>

#if !defined(BOOST_NO_CXX11_HDR_TUPLE)

#include <tuple>

namespace user
{

struct Y1
{
    int a;
    int b;
};

template<std::size_t I> int& get( Y1& v );
template<std::size_t I> int const& get( Y1 const& v );

template<> int& get<0>( Y1& v )
{
    return v.a;
}

template<> int const& get<0>( Y1 const& v )
{
    return v.a;
}

template<> int& get<1>( Y1& v )
{
    return v.b;
}

template<> int const& get<1>( Y1 const& v )
{
    return v.b;
}

struct Y2
{
    int a;
    int b;

    template<class T> friend
        typename boost::enable_if_<boost::is_same<T, Y2>::value, std::size_t>::type
        hash_value( T const& v )
    {
        std::size_t seed = 0;

        boost::hash_combine( seed, v.a );
        boost::hash_combine( seed, v.b );

        return seed;
    }
};

} // namespace user

namespace std
{

template<> struct tuple_size<user::Y1>: std::integral_constant<std::size_t, 2>
{
};

template<> struct tuple_size<user::Y2>: std::integral_constant<std::size_t, 2>
{
};

} // namespace std

namespace boost
{
namespace container_hash
{

    template<> struct is_tuple_like<user::Y2>: boost::false_type {};

} // namespace container_hash
} // namespace boost

#endif

template<class T> std::size_t hv( T const& t )
{
    return boost::hash<T>()( t );
}

int main()
{
    {
        std::pair<int, int> tp( 1, 2 );
        int const a[] = { 1, 2 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

#if !defined(BOOST_NO_CXX11_HDR_TUPLE)

    {
        std::tuple<> tp;

        BOOST_TEST_EQ( hv(tp), 0u );
    }

    {
        std::tuple<int> tp( 1 );
        int const a[] = { 1 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int> tp( 1, 2 );
        int const a[] = { 1, 2 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int, int> tp( 1, 2, 3 );
        int const a[] = { 1, 2, 3 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int, int, int> tp( 1, 2, 3, 4 );
        int const a[] = { 1, 2, 3, 4 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int, int, int, int> tp( 1, 2, 3, 4, 5 );
        int const a[] = { 1, 2, 3, 4, 5 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int, int, int, int, int> tp( 1, 2, 3, 4, 5, 6 );
        int const a[] = { 1, 2, 3, 4, 5, 6 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int, int, int, int, int, int> tp( 1, 2, 3, 4, 5, 6, 7 );
        int const a[] = { 1, 2, 3, 4, 5, 6, 7 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int, int, int, int, int, int, int> tp( 1, 2, 3, 4, 5, 6, 7, 8 );
        int const a[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

    {
        std::tuple<int, int, int, int, int, int, int, int, int> tp( 1, 2, 3, 4, 5, 6, 7, 8, 9 );
        int const a[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

#if !BOOST_WORKAROUND(BOOST_MSVC, <= 1800)

    {
        user::Y1 tp = { 1, 2 };
        int const a[] = { 1, 2 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

#endif

    {
        user::Y2 tp = { 1, 2 };
        int const a[] = { 1, 2 };

        BOOST_TEST_EQ( hv(tp), hv(a) );
    }

#endif

    return boost::report_errors();
}
