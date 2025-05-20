// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__clang__)
# pragma clang diagnostic ignored "-Wunused-private-field"
#endif

#include <boost/container_hash/hash.hpp>
#include <boost/describe/class.hpp>
#include <boost/describe/operators.hpp>
#include <boost/core/lightweight_test_trait.hpp>

#if !defined(BOOST_DESCRIBE_CXX14)

#include <boost/config/pragma_message.hpp>

BOOST_PRAGMA_MESSAGE( "Skipping test because BOOST_DESCRIBE_CXX14 is not defined" )

int main() {}

#else

struct X1
{
    int m;
    explicit X1( int m_ ): m( m_ ) {}
};

BOOST_DESCRIBE_STRUCT( X1, (), (m) )

struct X2
{
    int m;
    explicit X2( int m_ ): m( m_ ) {}
};

BOOST_DESCRIBE_STRUCT( X2, (), (m) )

struct X3
{
    int m;
    explicit X3( int m_ ): m( m_ ) {}
};

BOOST_DESCRIBE_STRUCT( X3, (), (m) )

class Y: public X1, protected X2, private X3
{
public:

    int m1;

protected:

    int m2;

private:

    int m3;

public:

    Y( int x1, int x2, int x3, int m1_, int m2_, int m3_ ):
      X1( x1 ), X2( x2 ), X3( x3 ), m1( m1_ ), m2( m2_ ), m3( m3_ ) {}

    BOOST_DESCRIBE_CLASS( Y, (X1, X2, X3), (m1), (m2), (m3) )
};

using boost::describe::operators::operator==;
using boost::describe::operators::operator!=;
using boost::describe::operators::operator<<;

int main()
{
    Y y1( 0, 0, 0, 0, 0, 0 );

    BOOST_TEST_EQ( y1, y1 );
    BOOST_TEST_EQ( boost::hash<Y>()(y1), boost::hash<Y>()(y1) );

    Y y2( 1, 0, 0, 0, 0, 0 );

    BOOST_TEST_NE( y1, y2 );
    BOOST_TEST_NE( boost::hash<Y>()(y1), boost::hash<Y>()(y2) );

    Y y3( 0, 1, 0, 0, 0, 0 );

    BOOST_TEST_NE( y1, y3 );
    BOOST_TEST_NE( boost::hash<Y>()(y1), boost::hash<Y>()(y3) );

    BOOST_TEST_NE( y2, y3 );
    BOOST_TEST_NE( boost::hash<Y>()(y2), boost::hash<Y>()(y3) );

    Y y4( 0, 0, 1, 0, 0, 0 );

    BOOST_TEST_NE( y1, y4 );
    BOOST_TEST_NE( boost::hash<Y>()(y1), boost::hash<Y>()(y4) );

    BOOST_TEST_NE( y2, y4 );
    BOOST_TEST_NE( boost::hash<Y>()(y2), boost::hash<Y>()(y4) );

    BOOST_TEST_NE( y3, y4 );
    BOOST_TEST_NE( boost::hash<Y>()(y3), boost::hash<Y>()(y4) );

    Y y5( 0, 0, 0, 1, 0, 0 );

    BOOST_TEST_NE( y1, y5 );
    BOOST_TEST_NE( boost::hash<Y>()(y1), boost::hash<Y>()(y5) );

    BOOST_TEST_NE( y2, y5 );
    BOOST_TEST_NE( boost::hash<Y>()(y2), boost::hash<Y>()(y5) );

    BOOST_TEST_NE( y3, y5 );
    BOOST_TEST_NE( boost::hash<Y>()(y3), boost::hash<Y>()(y5) );

    BOOST_TEST_NE( y4, y5 );
    BOOST_TEST_NE( boost::hash<Y>()(y4), boost::hash<Y>()(y5) );

    Y y6( 0, 0, 0, 0, 1, 0 );

    BOOST_TEST_NE( y1, y6 );
    BOOST_TEST_NE( boost::hash<Y>()(y1), boost::hash<Y>()(y6) );

    BOOST_TEST_NE( y2, y6 );
    BOOST_TEST_NE( boost::hash<Y>()(y2), boost::hash<Y>()(y6) );

    BOOST_TEST_NE( y3, y6 );
    BOOST_TEST_NE( boost::hash<Y>()(y3), boost::hash<Y>()(y6) );

    BOOST_TEST_NE( y4, y6 );
    BOOST_TEST_NE( boost::hash<Y>()(y4), boost::hash<Y>()(y6) );

    BOOST_TEST_NE( y5, y6 );
    BOOST_TEST_NE( boost::hash<Y>()(y5), boost::hash<Y>()(y6) );

    Y y7( 0, 0, 0, 0, 0, 1 );

    BOOST_TEST_NE( y1, y7 );
    BOOST_TEST_NE( boost::hash<Y>()(y1), boost::hash<Y>()(y7) );

    BOOST_TEST_NE( y2, y7 );
    BOOST_TEST_NE( boost::hash<Y>()(y2), boost::hash<Y>()(y7) );

    BOOST_TEST_NE( y3, y7 );
    BOOST_TEST_NE( boost::hash<Y>()(y3), boost::hash<Y>()(y7) );

    BOOST_TEST_NE( y4, y7 );
    BOOST_TEST_NE( boost::hash<Y>()(y4), boost::hash<Y>()(y7) );

    BOOST_TEST_NE( y5, y7 );
    BOOST_TEST_NE( boost::hash<Y>()(y5), boost::hash<Y>()(y7) );

    BOOST_TEST_NE( y6, y7 );
    BOOST_TEST_NE( boost::hash<Y>()(y6), boost::hash<Y>()(y7) );

    return boost::report_errors();
}

#endif
