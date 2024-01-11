// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(_MSC_VER) && _MSC_VER < 1910
# pragma warning( disable: 4800 ) // forcing value to bool 'true' or 'false'
#endif

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <utility>
#include <type_traits>

using namespace boost::system;

// Tricky mixed construction cases
// https://github.com/boostorg/system/issues/104
// https://brevzin.github.io//c++/2023/01/18/optional-construction/

template<class R1, class R2> void test()
{
    {
        R1 r1( make_error_code( errc::invalid_argument ) );
        R2 r2( r1 );

        BOOST_TEST( !r2.has_value() );
    }

    {
        R1 r1( 0 );
        R2 r2( r1 );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), false );
    }

    {
        R1 r1( 1 );
        R2 r2( r1 );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), true );
    }

    {
        R1 r1( make_error_code( errc::invalid_argument ) );
        R2 r2( std::move( r1 ) );

        BOOST_TEST( !r2.has_value() );
    }

    {
        R1 r1( 0 );
        R2 r2( std::move( r1 ) );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), false );
    }

    {
        R1 r1( 1 );
        R2 r2( std::move( r1 ) );

        BOOST_TEST( r2.has_value() ) && BOOST_TEST_EQ( r2.value(), true );
    }
}

struct X
{
};

int main()
{
    test< result<int>, result<bool> >();
    test< result<int> const, result<bool> >();

    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool>, result<X>&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool>, result<X> const&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool>, result<X>&&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool>, result<X> const&&>));

    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool const>, result<X>&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool const>, result<X> const&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool const>, result<X>&&>));
    BOOST_TEST_TRAIT_FALSE((std::is_constructible<result<bool const>, result<X> const&&>));

    return boost::report_errors();
}
