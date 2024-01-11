// Copyright 2023 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/mp11/integral.hpp>

#if !defined(BOOST_MP11_HAS_TEMPLATE_AUTO)

#pragma message("Test skipped because BOOST_MP11_HAS_TEMPLATE_AUTO is not defined")
int main() {}

#else

#include <boost/core/lightweight_test_trait.hpp>
#include <type_traits>
#include <cstddef>

int main()
{
    using boost::mp11::mp_value;

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<false>, std::integral_constant<bool, false>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<true>, std::integral_constant<bool, true>>));

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<'a'>, std::integral_constant<char, 'a'>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<'0'>, std::integral_constant<char, '0'>>));

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<0>, std::integral_constant<int, 0>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<-1>, std::integral_constant<int, -1>>));

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<0L>, std::integral_constant<long, 0>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<-1L>, std::integral_constant<long, -1>>));

    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<0u>, std::integral_constant<unsigned, 0>>));
    BOOST_TEST_TRAIT_TRUE((std::is_same<mp_value<1u>, std::integral_constant<unsigned, 1>>));

    return boost::report_errors();
}

#endif
