// Copyright 2023 Braden Ganetsky (braden.ganetsky@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
//
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test_trait.hpp>
#include <boost/mp11/algorithm.hpp>
#include <type_traits>
#include <tuple>
#include <utility>

int main()
{
    using boost::mp11::mp_list;
    using boost::mp11::mp_size_t;
    using boost::mp11::mp_slice_c;
    using boost::mp11::mp_slice;

    {
        using L1 = mp_list<int[], void, float&, char, bool*>;

        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 0, 0>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 0, 1>, mp_list<int[]>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 0, 2>, mp_list<int[], void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 0, 3>, mp_list<int[], void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 0, 4>, mp_list<int[], void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 0, 5>, mp_list<int[], void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 1, 1>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 1, 2>, mp_list<void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 1, 3>, mp_list<void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 1, 4>, mp_list<void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 1, 5>, mp_list<void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 2, 2>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 2, 3>, mp_list<float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 2, 4>, mp_list<float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 2, 5>, mp_list<float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 3, 3>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 3, 4>, mp_list<char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 3, 5>, mp_list<char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 4, 4>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 4, 5>, mp_list<bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L1, 5, 5>, mp_list<>>));

        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<0>, mp_size_t<0>>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<0>, mp_size_t<1>>, mp_list<int[]>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<0>, mp_size_t<2>>, mp_list<int[], void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<0>, mp_size_t<3>>, mp_list<int[], void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<0>, mp_size_t<4>>, mp_list<int[], void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<0>, mp_size_t<5>>, mp_list<int[], void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<1>, mp_size_t<1>>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<1>, mp_size_t<2>>, mp_list<void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<1>, mp_size_t<3>>, mp_list<void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<1>, mp_size_t<4>>, mp_list<void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<1>, mp_size_t<5>>, mp_list<void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<2>, mp_size_t<2>>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<2>, mp_size_t<3>>, mp_list<float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<2>, mp_size_t<4>>, mp_list<float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<2>, mp_size_t<5>>, mp_list<float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<3>, mp_size_t<3>>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<3>, mp_size_t<4>>, mp_list<char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<3>, mp_size_t<5>>, mp_list<char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<4>, mp_size_t<4>>, mp_list<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<4>, mp_size_t<5>>, mp_list<bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L1, mp_size_t<5>, mp_size_t<5>>, mp_list<>>));
    }

    {
        using L2 = std::tuple<int[], void, float&, char, bool*>;

        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 0, 0>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 0, 1>, std::tuple<int[]>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 0, 2>, std::tuple<int[], void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 0, 3>, std::tuple<int[], void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 0, 4>, std::tuple<int[], void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 0, 5>, std::tuple<int[], void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 1, 1>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 1, 2>, std::tuple<void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 1, 3>, std::tuple<void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 1, 4>, std::tuple<void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 1, 5>, std::tuple<void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 2, 2>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 2, 3>, std::tuple<float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 2, 4>, std::tuple<float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 2, 5>, std::tuple<float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 3, 3>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 3, 4>, std::tuple<char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 3, 5>, std::tuple<char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 4, 4>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 4, 5>, std::tuple<bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice_c<L2, 5, 5>, std::tuple<>>));

        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<0>, mp_size_t<0>>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<0>, mp_size_t<1>>, std::tuple<int[]>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<0>, mp_size_t<2>>, std::tuple<int[], void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<0>, mp_size_t<3>>, std::tuple<int[], void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<0>, mp_size_t<4>>, std::tuple<int[], void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<0>, mp_size_t<5>>, std::tuple<int[], void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<1>, mp_size_t<1>>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<1>, mp_size_t<2>>, std::tuple<void>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<1>, mp_size_t<3>>, std::tuple<void, float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<1>, mp_size_t<4>>, std::tuple<void, float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<1>, mp_size_t<5>>, std::tuple<void, float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<2>, mp_size_t<2>>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<2>, mp_size_t<3>>, std::tuple<float&>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<2>, mp_size_t<4>>, std::tuple<float&, char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<2>, mp_size_t<5>>, std::tuple<float&, char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<3>, mp_size_t<3>>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<3>, mp_size_t<4>>, std::tuple<char>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<3>, mp_size_t<5>>, std::tuple<char, bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<4>, mp_size_t<4>>, std::tuple<>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<4>, mp_size_t<5>>, std::tuple<bool*>>));
        BOOST_TEST_TRAIT_TRUE((std::is_same<mp_slice<L2, mp_size_t<5>, mp_size_t<5>>, std::tuple<>>));
    }

    return boost::report_errors();
}
