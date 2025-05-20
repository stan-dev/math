//  Copyright (C) 2011 Tim Blechmann
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <boost/lockfree/stack.hpp>

#define BOOST_TEST_MAIN
#ifdef BOOST_LOCKFREE_INCLUDE_TESTS
#    include <boost/test/included/unit_test.hpp>
#else
#    include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE( simple_stack_test )
{
    boost::lockfree::stack< long > stk( 128 );

    stk.push( 1 );
    stk.push( 2 );
    long out;
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.pop( out ) );
}

BOOST_AUTO_TEST_CASE( unsafe_stack_test )
{
    boost::lockfree::stack< long > stk( 128 );

    stk.unsynchronized_push( 1 );
    stk.unsynchronized_push( 2 );
    long out;
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.unsynchronized_pop( out ) );
}

BOOST_AUTO_TEST_CASE( ranged_push_test )
{
    boost::lockfree::stack< long > stk( 128 );

    long data[ 2 ] = { 1, 2 };

    BOOST_TEST_REQUIRE( stk.push( data, data + 2 ) == data + 2 );

    long out;
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.unsynchronized_pop( out ) );
}

BOOST_AUTO_TEST_CASE( span_push_test )
{
    boost::lockfree::stack< long > stk( 128 );

    long data[ 2 ] = { 1, 2 };

    BOOST_TEST_REQUIRE( stk.push( boost::span< const long >( data ) ) == 2 );

    long out;
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.unsynchronized_pop( out ) );
}


BOOST_AUTO_TEST_CASE( ranged_unsynchronized_push_test )
{
    boost::lockfree::stack< long > stk( 128 );

    long data[ 2 ] = { 1, 2 };

    BOOST_TEST_REQUIRE( stk.unsynchronized_push( data, data + 2 ) == data + 2 );

    long out;
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.unsynchronized_pop( out ) );
}

BOOST_AUTO_TEST_CASE( span_unsynchronized_push_test )
{
    boost::lockfree::stack< long > stk( 128 );

    long data[ 2 ] = { 1, 2 };

    BOOST_TEST_REQUIRE( stk.unsynchronized_push( boost::span< const long >( data ) ) == 2 );

    long out;
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.unsynchronized_pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.unsynchronized_pop( out ) );
}


BOOST_AUTO_TEST_CASE( fixed_size_stack_test )
{
    boost::lockfree::stack< long, boost::lockfree::capacity< 128 > > stk;

    stk.push( 1 );
    stk.push( 2 );
    long out;
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.pop( out ) );
    BOOST_TEST_REQUIRE( stk.empty() );
}

BOOST_AUTO_TEST_CASE( fixed_size_stack_test_exhausted )
{
    boost::lockfree::stack< long, boost::lockfree::capacity< 2 > > stk;

    stk.push( 1 );
    stk.push( 2 );
    BOOST_TEST_REQUIRE( !stk.push( 3 ) );
    long out;
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.pop( out ) );
    BOOST_TEST_REQUIRE( stk.empty() );
}

BOOST_AUTO_TEST_CASE( bounded_stack_test_exhausted )
{
    boost::lockfree::stack< long > stk( 2 );

    stk.bounded_push( 1 );
    stk.bounded_push( 2 );
    BOOST_TEST_REQUIRE( !stk.bounded_push( 3 ) );
    long out;
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 2 );
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( out == 1 );
    BOOST_TEST_REQUIRE( !stk.pop( out ) );
    BOOST_TEST_REQUIRE( stk.empty() );
}

BOOST_AUTO_TEST_CASE( stack_consume_one_test )
{
    boost::lockfree::stack< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );
    BOOST_TEST_REQUIRE( f.empty() );

    f.push( 1 );
    f.push( 2 );

    bool success1 = f.consume_one( []( int i ) {
        BOOST_TEST_REQUIRE( i == 2 );
    } );

    bool success2 = f.consume_one( []( int i ) mutable {
        BOOST_TEST_REQUIRE( i == 1 );
    } );

    BOOST_TEST_REQUIRE( success1 );
    BOOST_TEST_REQUIRE( success2 );

    BOOST_TEST_REQUIRE( f.empty() );
}

BOOST_AUTO_TEST_CASE( stack_consume_all_test )
{
    boost::lockfree::stack< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );
    BOOST_TEST_REQUIRE( f.empty() );

    f.push( 1 );
    f.push( 2 );

    size_t consumed = f.consume_all( []( int i ) {} );

    BOOST_TEST_REQUIRE( consumed == 2u );

    BOOST_TEST_REQUIRE( f.empty() );
}

BOOST_AUTO_TEST_CASE( stack_consume_all_atomic_test )
{
    boost::lockfree::stack< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );
    BOOST_TEST_REQUIRE( f.empty() );

    f.push( 1 );
    f.push( 2 );
    f.push( 3 );

    size_t consumed = f.consume_all_atomic( []( int i ) {} );

    BOOST_TEST_REQUIRE( consumed == 3u );

    BOOST_TEST_REQUIRE( f.empty() );
}


BOOST_AUTO_TEST_CASE( stack_consume_all_atomic_reversed_test )
{
    boost::lockfree::stack< int > f( 64 );

    BOOST_TEST_WARN( f.is_lock_free() );
    BOOST_TEST_REQUIRE( f.empty() );

    f.push( 1 );
    f.push( 2 );
    f.push( 3 );

    size_t consumed = f.consume_all_atomic_reversed( []( int i ) {} );

    BOOST_TEST_REQUIRE( consumed == 3u );

    BOOST_TEST_REQUIRE( f.empty() );
}


BOOST_AUTO_TEST_CASE( reserve_test )
{
    typedef boost::lockfree::stack< void* > memory_stack;

    memory_stack ms( 1 );
    ms.reserve( 1 );
    ms.reserve_unsafe( 1 );
}

BOOST_AUTO_TEST_CASE( stack_with_allocator )
{
    using allocator_type = std::allocator< char >;

    using stack_t = boost::lockfree::stack< char, boost::lockfree::allocator< allocator_type > >;
    using stack_with_capacity_t
        = boost::lockfree::stack< char, boost::lockfree::allocator< allocator_type >, boost::lockfree::capacity< 16 > >;

    auto allocator = stack_t::allocator {};

    {
        stack_with_capacity_t stack_with_allocator {
            allocator,
        };
        stack_t stack_with_size_and_allocator {
            5,
            allocator,
        };
    }
    {
        stack_with_capacity_t stack_with_allocator {
            allocator_type {},
        };
        stack_t stack_with_size_and_allocator {
            5,
            allocator_type {},
        };
    }
}

BOOST_AUTO_TEST_CASE( move_semantics )
{
    boost::lockfree::stack< std::unique_ptr< int >, boost::lockfree::capacity< 128 > > stk;

    stk.push( std::make_unique< int >( 0 ) );
    stk.push( std::make_unique< int >( 1 ) );

    auto two = std::make_unique< int >( 2 );
    stk.push( std::move( two ) );

    std::unique_ptr< int > out;
    BOOST_TEST_REQUIRE( stk.pop( out ) );
    BOOST_TEST_REQUIRE( *out == 2 );

    stk.consume_one( []( std::unique_ptr< int > one ) {
        BOOST_TEST_REQUIRE( *one == 1 );
    } );

    stk.consume_all( []( std::unique_ptr< int > ) {} );
}

#if !defined( BOOST_NO_CXX17_HDR_OPTIONAL )

BOOST_AUTO_TEST_CASE( queue_uses_optional )
{
    boost::lockfree::stack< int > stk( 5 );

    bool pop_to_nullopt = stk.pop( boost::lockfree::uses_optional ) == std::nullopt;
    BOOST_TEST_REQUIRE( pop_to_nullopt );

    stk.push( 53 );
    bool pop_to_optional = stk.pop( boost::lockfree::uses_optional ) == 53;
    BOOST_TEST_REQUIRE( pop_to_optional );
}

#endif
