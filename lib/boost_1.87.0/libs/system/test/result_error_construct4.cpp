// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system/result.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::system;

// Eigen::Matrix4d has an explicit templated constructor
// https://github.com/boostorg/system/issues/103
// https://github.com/boostorg/json/issues/843

struct X
{
    X() {}
    template<class T> explicit X( T const& ) {}
};

int main()
{
    {
        auto ec = make_error_code( errc::invalid_argument );

        result<X> r = ec;

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), ec );
    }

    {
        auto ec = make_error_code( std::errc::invalid_argument );

        result<X> r = ec;

        BOOST_TEST( !r.has_value() );
        BOOST_TEST( r.has_error() );

        BOOST_TEST_EQ( r.error(), ec );
    }

    return boost::report_errors();
}
