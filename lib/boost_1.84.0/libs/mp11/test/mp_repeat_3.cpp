// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/detail/config.hpp>

#if BOOST_MP11_MSVC
# pragma warning( disable: 4503 ) // decorated name length exceeded
#endif

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/mp11/function.hpp>
#include <boost/core/lightweight_test.hpp>

int main()
{
    using boost::mp11::mp_repeat_c;
    using boost::mp11::mp_list;
    using boost::mp11::mp_int;
    using boost::mp11::mp_plus;
    using boost::mp11::mp_apply;

    int const N = 529;

    {
        using L1 = mp_repeat_c< mp_list<mp_int<1>, mp_int<2>>, N >;
        using R1 = mp_apply<mp_plus, L1>;

        BOOST_TEST_EQ( R1::value, N * 3 );
    }

#if defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

    using boost::mp11::mp_list_v;

    {
        using L1 = mp_repeat_c< mp_list_v<1, 2>, N >;
        using R1 = mp_apply<mp_plus, L1>;

        BOOST_TEST_EQ( R1::value, N * 3 );
    }

#endif

    return boost::report_errors();
}
