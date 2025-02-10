//  Copyright (C) 2011 Tim Blechmann
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <boost/lockfree/lockfree_forward.hpp>

#include <boost/lockfree/queue.hpp>
#include <boost/thread.hpp>

#define BOOST_TEST_MAIN
#ifdef BOOST_LOCKFREE_INCLUDE_TESTS
#    include <boost/test/included/unit_test.hpp>
#else
#    include <boost/test/unit_test.hpp>
#endif

#include <memory>


using namespace boost;
using namespace boost::lockfree;
using namespace std;

BOOST_AUTO_TEST_CASE( simple_queue_test )
{
    queue< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );

    BOOST_TEST_REQUIRE( f.empty() );
    f.push( 1 );
    f.push( 2 );

    int i1( 0 ), i2( 0 );

    BOOST_TEST_REQUIRE( f.pop( i1 ) );
    BOOST_TEST_REQUIRE( i1 == 1 );

    BOOST_TEST_REQUIRE( f.pop( i2 ) );
    BOOST_TEST_REQUIRE( i2 == 2 );
    BOOST_TEST_REQUIRE( f.empty() );
}

BOOST_AUTO_TEST_CASE( simple_queue_test_capacity )
{
    queue< int, capacity< 64 > > f;

    BOOST_TEST_WARN( f.is_lock_free() );

    BOOST_TEST_REQUIRE( f.empty() );
    f.push( 1 );
    f.push( 2 );

    int i1( 0 ), i2( 0 );

    BOOST_TEST_REQUIRE( f.pop( i1 ) );
    BOOST_TEST_REQUIRE( i1 == 1 );

    BOOST_TEST_REQUIRE( f.pop( i2 ) );
    BOOST_TEST_REQUIRE( i2 == 2 );
    BOOST_TEST_REQUIRE( f.empty() );
}


BOOST_AUTO_TEST_CASE( unsafe_queue_test )
{
    queue< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );
    BOOST_TEST_REQUIRE( f.empty() );

    int i1( 0 ), i2( 0 );

    f.unsynchronized_push( 1 );
    f.unsynchronized_push( 2 );

    BOOST_TEST_REQUIRE( f.unsynchronized_pop( i1 ) );
    BOOST_TEST_REQUIRE( i1 == 1 );

    BOOST_TEST_REQUIRE( f.unsynchronized_pop( i2 ) );
    BOOST_TEST_REQUIRE( i2 == 2 );
    BOOST_TEST_REQUIRE( f.empty() );
}


BOOST_AUTO_TEST_CASE( queue_consume_one_test )
{
    queue< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );
    BOOST_TEST_REQUIRE( f.empty() );

    f.push( 1 );
    f.push( 2 );

    bool success1 = f.consume_one( []( int i ) {
        BOOST_TEST_REQUIRE( i == 1 );
    } );

    bool success2 = f.consume_one( []( int i ) mutable {
        BOOST_TEST_REQUIRE( i == 2 );
    } );

    BOOST_TEST_REQUIRE( success1 );
    BOOST_TEST_REQUIRE( success2 );

    BOOST_TEST_REQUIRE( f.empty() );
}

BOOST_AUTO_TEST_CASE( queue_consume_all_test )
{
    queue< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );
    BOOST_TEST_REQUIRE( f.empty() );

    f.push( 1 );
    f.push( 2 );

    size_t consumed = f.consume_all( []( int i ) {} );

    BOOST_TEST_REQUIRE( consumed == 2u );

    BOOST_TEST_REQUIRE( f.empty() );
}


BOOST_AUTO_TEST_CASE( queue_convert_pop_test )
{
    queue< int* > f( 128 );
    BOOST_TEST_REQUIRE( f.empty() );
    f.push( new int( 1 ) );
    f.push( new int( 2 ) );
    f.push( new int( 3 ) );
    f.push( new int( 4 ) );

    {
        int* i1;

        BOOST_TEST_REQUIRE( f.pop( i1 ) );
        BOOST_TEST_REQUIRE( *i1 == 1 );
        delete i1;
    }


    {
        boost::shared_ptr< int > i2;
        BOOST_TEST_REQUIRE( f.pop( i2 ) );
        BOOST_TEST_REQUIRE( *i2 == 2 );
    }

    {
        unique_ptr< int > i3;
        BOOST_TEST_REQUIRE( f.pop( i3 ) );

        BOOST_TEST_REQUIRE( *i3 == 3 );
    }

    {
        std::shared_ptr< int > i4;
        BOOST_TEST_REQUIRE( f.pop( i4 ) );

        BOOST_TEST_REQUIRE( *i4 == 4 );
    }


    BOOST_TEST_REQUIRE( f.empty() );
}

BOOST_AUTO_TEST_CASE( reserve_test )
{
    typedef boost::lockfree::queue< void* > memory_queue;

    memory_queue ms( 1 );
    ms.reserve( 1 );
    ms.reserve_unsafe( 1 );
}

BOOST_AUTO_TEST_CASE( queue_with_allocator )
{
    using allocator_type = std::allocator< char >;

    using queue_t = boost::lockfree::queue< char, boost::lockfree::allocator< allocator_type > >;
    using queue_with_capacity_t
        = boost::lockfree::queue< char, boost::lockfree::allocator< allocator_type >, boost::lockfree::capacity< 16 > >;

    auto allocator = queue_t::allocator {};

    {
        queue_with_capacity_t q_with_allocator {
            allocator,
        };
        queue_t q_with_size_and_allocator {
            5,
            allocator,
        };
    }
    {
        queue_with_capacity_t q_with_allocator {
            allocator_type {},
        };
        queue_t q_with_size_and_allocator {
            5,
            allocator_type {},
        };
    }
}

BOOST_AUTO_TEST_CASE( move_semantics )
{
    boost::lockfree::queue< int, boost::lockfree::capacity< 128 > > stk;

    stk.push( 0 );
    stk.push( 1 );

    auto two = 2;
    stk.push( std::move( two ) );

    int out;
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 0 );

    stk.consume_one( []( int one ) {
        BOOST_TEST_REQUIRE( one == 1 );
    } );

    stk.consume_all( []( int ) {} );
}

#if !defined( BOOST_NO_CXX17_HDR_OPTIONAL )

BOOST_AUTO_TEST_CASE( queue_uses_optional )
{
    boost::lockfree::queue< int > stk( 5 );

    bool pop_to_nullopt = stk.pop( boost::lockfree::uses_optional ) == std::nullopt;
    BOOST_TEST_REQUIRE( pop_to_nullopt );

    stk.push( 53 );
    bool pop_to_optional = stk.pop( boost::lockfree::uses_optional ) == 53;
    BOOST_TEST_REQUIRE( pop_to_optional );
}

#endif
