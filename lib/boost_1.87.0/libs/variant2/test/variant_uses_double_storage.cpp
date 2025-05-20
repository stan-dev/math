// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/variant2/variant.hpp>
#include <type_traits>

using namespace boost::variant2;

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)

struct X1
{
};

STATIC_ASSERT( std::is_nothrow_move_constructible<X1>::value );
STATIC_ASSERT( std::is_trivially_destructible<X1>::value );

struct X2
{
    ~X2() {}
};

STATIC_ASSERT( std::is_nothrow_move_constructible<X2>::value );
STATIC_ASSERT( !std::is_trivially_destructible<X2>::value );

struct X3
{
    X3( X3&& ) {}
};

STATIC_ASSERT( !std::is_nothrow_move_constructible<X3>::value );
STATIC_ASSERT( std::is_trivially_destructible<X3>::value );

struct X4
{
    ~X4() {}
    X4( X4&& ) {}
};

STATIC_ASSERT( !std::is_nothrow_move_constructible<X4>::value );
STATIC_ASSERT( !std::is_trivially_destructible<X4>::value );

//

STATIC_ASSERT( !variant<int, float>::uses_double_storage() );
STATIC_ASSERT( !variant<int, float, X1>::uses_double_storage() );
STATIC_ASSERT( !variant<int, float, X2>::uses_double_storage() );
STATIC_ASSERT( variant<int, float, X3>::uses_double_storage() );
STATIC_ASSERT( variant<int, float, X4>::uses_double_storage() );
