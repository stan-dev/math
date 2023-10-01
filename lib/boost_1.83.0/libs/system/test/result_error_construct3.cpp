// Copyright 2017, 2021, 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>
#include <cerrno>

using namespace boost::system;

struct X
{
    static int instances;

    int v_;

    X(): v_() { ++instances; }

    explicit X( int v ): v_( v ) { ++instances; }

    X( int v1, int v2 ): v_( v1+v2 ) { ++instances; }
    X( int v1, int v2, int v3 ): v_( v1+v2+v3 ) { ++instances; }

    X( X const& r ): v_( r.v_ ) { ++instances; }

    X& operator=( X const& ) = delete;

    ~X() { --instances; }
};

int X::instances = 0;

int main()
{
    {
        auto ec = make_error_code( errc::invalid_argument );

        using R = result<int>;
        R r( R::in_place_error, ec );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), ec );
    }

    {
        using R = result<int>;
        R r( R::in_place_error, EINVAL, generic_category() );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), error_code( EINVAL, generic_category() ) );
    }

    {
        auto ec = make_error_code( errc::invalid_argument );

        using R = result<error_code>;
        R r( R::in_place_error, ec );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), ec );
    }

    {
        using R = result<error_code>;
        R r( R::in_place_error, EINVAL, generic_category() );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), error_code( EINVAL, generic_category() ) );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<std::string, X>;
        R r( R::in_place_error, 1 );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error().v_, 1 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<int, X>;
        R r( R::in_place_error, 1, 2 );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error().v_, 1+2 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<int, X>;
        R r( R::in_place_error, 1, 2, 3 );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error().v_, 1+2+3 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        using R = result<X, X>;
        R r( R::in_place_error, 1 );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error().v_, 1 );

        BOOST_TEST_EQ( X::instances, 1 );
    }

    BOOST_TEST_EQ( X::instances, 0 );

    {
        auto ec = make_error_code( errc::invalid_argument );

        using R = result<void>;
        R r( R::in_place_error, ec );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), ec );
    }

    {
        using R = result<void>;
        R r( R::in_place_error, EINVAL, generic_category() );

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), error_code( EINVAL, generic_category() ) );
    }

    return boost::report_errors();
}
