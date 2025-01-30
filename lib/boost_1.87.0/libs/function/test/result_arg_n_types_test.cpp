// Copyright 2024 Peter Dimov
// Use, modification and distribution is subject to
// the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/function.hpp>
#include <boost/core/lightweight_test_trait.hpp>

struct R {};

struct A1 {};
struct A2 {};
struct A3 {};
struct A4 {};
struct A5 {};
struct A6 {};
struct A7 {};
struct A8 {};
struct A9 {};
struct A10 {};

int main()
{
    {
        typedef boost::function<R()> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
    }

    {
        typedef boost::function<R(A1)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
    }

    {
        typedef boost::function<R(A1, A2)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
    }

    {
        typedef boost::function<R(A1, A2, A3)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
    }

    {
        typedef boost::function<R(A1, A2, A3, A4)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
        BOOST_TEST_TRAIT_SAME(F::arg4_type, A4);
    }

    {
        typedef boost::function<R(A1, A2, A3, A4, A5)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
        BOOST_TEST_TRAIT_SAME(F::arg4_type, A4);
        BOOST_TEST_TRAIT_SAME(F::arg5_type, A5);
    }

    {
        typedef boost::function<R(A1, A2, A3, A4, A5, A6)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
        BOOST_TEST_TRAIT_SAME(F::arg4_type, A4);
        BOOST_TEST_TRAIT_SAME(F::arg5_type, A5);
        BOOST_TEST_TRAIT_SAME(F::arg6_type, A6);
    }

    {
        typedef boost::function<R(A1, A2, A3, A4, A5, A6, A7)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
        BOOST_TEST_TRAIT_SAME(F::arg4_type, A4);
        BOOST_TEST_TRAIT_SAME(F::arg5_type, A5);
        BOOST_TEST_TRAIT_SAME(F::arg6_type, A6);
        BOOST_TEST_TRAIT_SAME(F::arg7_type, A7);
    }

    {
        typedef boost::function<R(A1, A2, A3, A4, A5, A6, A7, A8)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
        BOOST_TEST_TRAIT_SAME(F::arg4_type, A4);
        BOOST_TEST_TRAIT_SAME(F::arg5_type, A5);
        BOOST_TEST_TRAIT_SAME(F::arg6_type, A6);
        BOOST_TEST_TRAIT_SAME(F::arg7_type, A7);
        BOOST_TEST_TRAIT_SAME(F::arg8_type, A8);
    }

    {
        typedef boost::function<R(A1, A2, A3, A4, A5, A6, A7, A8, A9)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
        BOOST_TEST_TRAIT_SAME(F::arg4_type, A4);
        BOOST_TEST_TRAIT_SAME(F::arg5_type, A5);
        BOOST_TEST_TRAIT_SAME(F::arg6_type, A6);
        BOOST_TEST_TRAIT_SAME(F::arg7_type, A7);
        BOOST_TEST_TRAIT_SAME(F::arg8_type, A8);
        BOOST_TEST_TRAIT_SAME(F::arg9_type, A9);
    }

    {
        typedef boost::function<R(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10)> F;

        BOOST_TEST_TRAIT_SAME(F::result_type, R);
        BOOST_TEST_TRAIT_SAME(F::arg1_type, A1);
        BOOST_TEST_TRAIT_SAME(F::arg2_type, A2);
        BOOST_TEST_TRAIT_SAME(F::arg3_type, A3);
        BOOST_TEST_TRAIT_SAME(F::arg4_type, A4);
        BOOST_TEST_TRAIT_SAME(F::arg5_type, A5);
        BOOST_TEST_TRAIT_SAME(F::arg6_type, A6);
        BOOST_TEST_TRAIT_SAME(F::arg7_type, A7);
        BOOST_TEST_TRAIT_SAME(F::arg8_type, A8);
        BOOST_TEST_TRAIT_SAME(F::arg9_type, A9);
        BOOST_TEST_TRAIT_SAME(F::arg10_type, A10);
    }

    return boost::report_errors();
}
