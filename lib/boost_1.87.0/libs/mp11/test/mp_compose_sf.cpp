// Copyright 2023 Jody Hagins
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11.hpp>
using namespace boost::mp11;

// A simple metafunction that takes the first three items from a list
template <typename T>
using take3 = mp_take<T, mp_int<3>>;

// Directly call take3
static_assert(mp_same<
    take3<mp_iota_c<9>>,
    mp_list<mp_size_t<0>, mp_size_t<1>, mp_size_t<2>>>::value, "");

// Invoke a quoted take3 (mp_quote<take3>)
static_assert(mp_same<
    mp_invoke_q<mp_quote<take3>, mp_iota_c<9>>,
    mp_list<mp_size_t<0>, mp_size_t<1>, mp_size_t<2>>>::value, "");

#if !BOOST_MP11_WORKAROUND( BOOST_MP11_MSVC, < 1900 )

// Invoke a composed take3 (mp_compose<take3>)
static_assert(mp_same<
    mp_invoke_q<mp_compose<take3>, mp_iota_c<9>>,
    mp_list<mp_size_t<0>, mp_size_t<1>, mp_size_t<2>>>::value, "");

#endif

// If we apply the following iteration, we expect to get all triplets of the list.
using expected = mp_list<
    mp_list<mp_size_t<0>, mp_size_t<1>, mp_size_t<2>>,
    mp_list<mp_size_t<3>, mp_size_t<4>, mp_size_t<5>>,
    mp_list<mp_size_t<6>, mp_size_t<7>, mp_size_t<8>>>;

// First using mp_quote - this works as expected
static_assert(
    mp_same<
        mp_iterate_q<
            mp_iota_c<9>,
            mp_quote<take3>,
            mp_bind_back<mp_drop, mp_int<3>>>,
        expected>::value, "");

#if !BOOST_MP11_WORKAROUND( BOOST_MP11_MSVC, < 1900 )

// Next, simply replace mp_quote with mp_compose.  We expect to get the same answer,
// but we get a hard compiler error instead.
static_assert(
    mp_same<
        mp_iterate_q<
            mp_iota_c<9>,
            mp_compose<take3>,
            mp_bind_back<mp_drop, mp_int<3>>>,
        expected>::value, "");

#endif
