// Copyright 2015, 2017, 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/mp11/detail/config.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <type_traits>

using boost::mp11::mp_apply_q;
using boost::mp11::mp_list;
using boost::mp11::mp_quote;
using boost::mp11::mp_quote_trait;
using boost::mp11::mp_valid;
using boost::mp11::mp_identity;
using boost::mp11::mp_identity_t;

int main()
{
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, void>));

#if !BOOST_MP11_WORKAROUND( BOOST_MP11_MSVC, < 1920 )

    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, void, void>));

    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, void, mp_list<>>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, void, mp_list<void>>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, void, mp_list<void, void>>));

#endif

    using Qi = mp_quote<mp_identity_t>;

    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, Qi, mp_list<>>));
    BOOST_TEST_TRAIT_TRUE((mp_valid<mp_apply_q, Qi, mp_list<void>>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, Qi, mp_list<void, void>>));

    using Qt = mp_quote_trait<mp_identity>;

#if !BOOST_MP11_WORKAROUND( BOOST_MP11_MSVC, < 1900 )
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, Qt, mp_list<>>));
#endif
    BOOST_TEST_TRAIT_TRUE((mp_valid<mp_apply_q, Qt, mp_list<void>>));
    BOOST_TEST_TRAIT_FALSE((mp_valid<mp_apply_q, Qt, mp_list<void, void>>));

    return boost::report_errors();
}
