// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/bind.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <utility>

struct X1 {};
struct X2 {};

struct Test
{
    void operator()( boost::mp11::mp_size_t<0> )
    {
    }

    template<class N> void operator()( N )
    {
        using boost::mp11::mp_repeat;
        using boost::mp11::mp_list;
        using boost::mp11::mp_apply_q;
        using boost::mp11::mp_bind_front;
        using boost::mp11::mp_transform_q;
        using boost::mp11::mp_quote;

        using L1 = mp_repeat< mp_list<std::pair<X1, X2>>, N >; // mp_list<pair<X1, X2>, pair<X1, X2>, ...>
        using L2 = mp_apply_q< mp_bind_front<mp_transform_q, mp_quote<mp_list>>, L1>; // pair<mp_list<X1, X1, ...>, mp_list<X2, X2, ...>>

        using R1 = mp_repeat<mp_list<X1>, N>; // mp_list<X1, X1, ...>
        using R2 = mp_repeat<mp_list<X2>, N>; // mp_list<X2, X2, ...>

        BOOST_TEST_TRAIT_SAME(L2, std::pair<R1, R2>);
    }
};

int main()
{
#if BOOST_MP11_WORKAROUND( BOOST_MP11_MSVC, < 1920 )

    using boost::mp11::mp_size_t;

    Test()( mp_size_t<1>() );
    Test()( mp_size_t<2>() );
    Test()( mp_size_t<3>() );
    Test()( mp_size_t<4>() );
    Test()( mp_size_t<5>() );
    Test()( mp_size_t<6>() );
    Test()( mp_size_t<7>() );
    Test()( mp_size_t<9>() );
    Test()( mp_size_t<10>() );
    Test()( mp_size_t<11>() );

#else

    using boost::mp11::mp_for_each;
    using boost::mp11::mp_iota_c;

    mp_for_each< mp_iota_c<32> >( Test() );

#endif

    return boost::report_errors();
}
