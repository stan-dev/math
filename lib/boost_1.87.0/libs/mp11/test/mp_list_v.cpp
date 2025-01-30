// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/list.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

#include <boost/core/lightweight_test.hpp>

template<class L> bool is_value_list( L const& ) { return false; }
template<template<auto...> class L, auto... A> bool is_value_list( L<A...> const& ) { return true; }

int main()
{
    using boost::mp11::mp_list_v;
    using boost::mp11::mp_list;

    BOOST_TEST_NOT(is_value_list( mp_list<>{} ));

    BOOST_TEST(is_value_list( mp_list_v<>{} ));
    BOOST_TEST(is_value_list( mp_list_v<false>{} ));
    BOOST_TEST(is_value_list( mp_list_v<false, 0>{} ));
    BOOST_TEST(is_value_list( mp_list_v<false, 0, std::size_t(0)>{} ));

    return boost::report_errors();
}

#endif
