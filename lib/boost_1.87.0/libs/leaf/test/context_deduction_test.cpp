// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/config.hpp>
#   include <boost/leaf/context.hpp>
#   include <boost/leaf/pred.hpp>
#   include <boost/leaf/diagnostics.hpp>
#endif

#include "_test_ec.hpp"

namespace leaf = boost::leaf;

template <class... T>
struct unwrap_tuple;

template <>
struct unwrap_tuple<std::tuple<>>
{
    using type = std::tuple<>;
};

template <template <class> class S, class... E>
struct unwrap_tuple<std::tuple<S<E>...>>
{
    using type = std::tuple<E...>;
};

template <class... H>
typename unwrap_tuple<typename std::decay<decltype(std::declval<typename leaf::context_type_from_handlers<H...>>().tup())>::type>::type * expd( H && ... )
{
    return 0;
}

template <class... T, class U>
void test( U * )
{
    static_assert(std::is_same<std::tuple<T...>, U>::value,"context deduction");
}

template <int> struct info { int value; };

enum class my_error_code
{
    ok,
    error1,
    error2,
    error3
};

namespace std
{
    template <> struct is_error_code_enum<my_error_code>: std::true_type { };
}

enum class my_error_condition
{
    ok,
    cond1
};

namespace std
{
    template <> struct is_error_condition_enum<my_error_condition>: std::true_type { };
}

void not_called_on_purpose()
{
    test<>( expd([]( leaf::error_info const & ){ }) );

    test<info<1>>( expd([]( info<1> ){ }) );
    test<info<1>>( expd([]( info<1> const ){ }) );
    test<info<1>>( expd([]( info<1> const & ){ }) );
    test<info<1>>( expd([]( info<1>, leaf::error_info const & ){ }) );

    test<info<1>, info<2>>( expd([]( info<1> ){ }, []( info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1> ){ }, []( info<1>, info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, info<2> ){ }, []( info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, info<2> ){ }, []( info<1>, info<2> ){ }) );

    test<info<1>, info<2>>( expd([]( leaf::error_info const &, info<1> ){ }, []( info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( leaf::error_info const &, info<1> ){ }, []( info<1>, info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( leaf::error_info const &, info<1>, info<2> ){ }, []( info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( leaf::error_info const &, info<1>, info<2> ){ }, []( info<1>, info<2> ){ }) );

    test<info<1>, info<2>>( expd([]( info<1>, leaf::error_info const & ){ }, []( info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, leaf::error_info const & ){ }, []( info<1>, info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, leaf::error_info const &, info<2> ){ }, []( info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, leaf::error_info const &, info<2> ){ }, []( info<1>, info<2> ){ }) );

    test<info<1>, info<2>>( expd([]( info<1> ){ }, []( leaf::error_info const &, info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1> ){ }, []( leaf::error_info const &, info<1>, info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, info<2> ){ }, []( leaf::error_info const &, info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, info<2> ){ }, []( leaf::error_info const &, info<1>, info<2> ){ }) );

    test<info<1>, info<2>>( expd([]( info<1> ){ }, []( info<2>, leaf::error_info const & ){ }) );
    test<info<1>, info<2>>( expd([]( info<1> ){ }, []( info<1>, leaf::error_info const &, info<2> ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, info<2> ){ }, []( info<2>, leaf::error_info const & ){ }) );
    test<info<1>, info<2>>( expd([]( info<1>, info<2> ){ }, []( info<1>, leaf::error_info const &, info<2> ){ }) );

    test<info<1>, info<2>, info<3>>( expd([]( info<1> ){ }, []( info<2> ){ }, []( info<3> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1> ){ }, []( info<1>, info<2> ){ }, []( info<1>, info<3> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1> ){ }, []( info<1>, info<2> ){ }, []( info<1>, info<3> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1>, info<2> ){ }, []( info<2> ){ }, []( info<3> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1>, info<2> ){ }, []( info<2> ){ }, []( info<3> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1> ){ }, []( info<2> ){ }, []( info<3>, info<2> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1> ){ }, []( info<2> ){ }, []( info<3>, info<2> ){ }) );
    test<info<1>, info<3>, info<2>>( expd([]( info<1>, info<3> ){ }, []( info<2> ){ }, []( info<3> ){ }) );
    test<info<1>, info<3>, info<2>>( expd([]( info<1>, info<3> ){ }, []( info<2> ){ }, []( info<3> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1> ){ }, []( info<2>, info<3> ){ }, []( info<3> ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1> ){ }, []( info<2>, info<3> ){ }, []( info<3> ){ }) );

    test<my_error_code>( expd([]( leaf::match<my_error_code, my_error_code::error1> ){ }) );
#if BOOST_LEAF_CFG_STD_SYSTEM_ERROR
    test<std::error_code>( expd([]( leaf::match<leaf::condition<cond_x>, cond_x::x00> ){ }) );
#endif

    test<info<1>>( expd([]( leaf::match_value<info<1>,42> ){ }) );

#if BOOST_LEAF_CFG_STD_SYSTEM_ERROR
    test<std::error_code>( expd([]( leaf::match<leaf::condition<my_error_condition>, my_error_condition::cond1> ){ }) );
#if __cplusplus >= 201703L
    test<std::error_code>( expd([]( leaf::match<std::error_code, my_error_code::error1> ){ }) );
    test<std::error_code>( expd([]( leaf::match<std::error_code, my_error_condition::cond1> ){ }) );
#endif
#endif

#ifndef BOOST_LEAF_NO_EXCEPTIONS
    test<info<1>>( expd([]( leaf::catch_<std::exception>, info<1> ){ }) );
#endif

    test<info<1>, info<2>, info<3>>( expd([]( info<1> const *, info<2> ){ }, []( info<1>, info<3> const * ){ }) );
    test<info<1>, info<2>, info<3>>( expd([]( info<1> const, info<2> ){ }, []( info<1> const *, info<3> ){ }) );

    test<info<1>, info<2>, leaf::e_source_location>( expd([]( info<1>, info<2>, leaf::diagnostic_info const & ){ }, []( info<1>, info<2> ){ }) );
#if BOOST_LEAF_CFG_DIAGNOSTICS && BOOST_LEAF_CFG_CAPTURE
    test<info<1>, info<2>, leaf::e_source_location, leaf::detail::dynamic_allocator>( expd([]( info<1>, info<2>, leaf::diagnostic_details const & ){ }, []( info<1>, info<2> ){ }) );
#else
    test<info<1>, info<2>, leaf::e_source_location>( expd([]( info<1>, info<2>, leaf::diagnostic_details const & ){ }, []( info<1>, info<2> ){ }) );
#endif

    {
        auto error_handlers = std::make_tuple(
            []( info<2> ) { },
            []( info<3> ) { } );
        test<info<1>, info<2>, info<3>, info<4>, info<5>>( expd(
            []( info<1> ) { },
            error_handlers,
            []( info<4>, info<5> ) { } ) );
    }

    {
        auto error_handlers1 = std::make_tuple(
            []( info<4> ) { },
            []( info<5> ) { } );
        auto error_handlers2 = std::make_tuple(
            []( info<2> ) { },
            []( info<3> ) { },
            error_handlers1 );
        test<info<1>, info<2>, info<3>, info<4>, info<5>, info<6>>( expd(
            []( info<1> ) { },
            error_handlers2,
            []( info<4>, info<5> ) { },
            []( info<6> ) { } ) );
    }
}

int main()
{
    return 0;
}
