//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil/extension/numeric/affine.hpp>
#include <boost/gil/extension/numeric/resample.hpp>
#include <boost/gil/extension/numeric/sampler.hpp>
#include <boost/gil.hpp>

#define BOOST_TEST_MODULE test_ext_numeric_matrix3x2
#include "unit_test.hpp"

namespace gil = boost::gil;

namespace {
constexpr double HALF_PI = 1.57079632679489661923;
}

BOOST_AUTO_TEST_CASE(matrix3x2_default_constructor)
{
    gil::matrix3x2<int> m1;
    BOOST_TEST(m1.a == 1);
    BOOST_TEST(m1.b == 0);
    BOOST_TEST(m1.c == 0);
    BOOST_TEST(m1.d == 1);
    BOOST_TEST(m1.e == 0);
    BOOST_TEST(m1.f == 0);
}

BOOST_AUTO_TEST_CASE(matrix3x2_parameterized_constructor)
{
    gil::matrix3x2<int> m1(1, 2, 3, 4, 5, 6);
    BOOST_TEST(m1.a == 1);
    BOOST_TEST(m1.b == 2);
    BOOST_TEST(m1.c == 3);
    BOOST_TEST(m1.d == 4);
    BOOST_TEST(m1.e == 5);
    BOOST_TEST(m1.f == 6);
}

BOOST_AUTO_TEST_CASE(matrix3x2_copy_constructor)
{
    gil::matrix3x2<int> m1(1, 2, 3, 4, 5, 6);
    gil::matrix3x2<int> m2(m1);
    BOOST_TEST(m2.a == 1);
    BOOST_TEST(m2.b == 2);
    BOOST_TEST(m2.c == 3);
    BOOST_TEST(m2.d == 4);
    BOOST_TEST(m2.e == 5);
    BOOST_TEST(m2.f == 6);
}

BOOST_AUTO_TEST_CASE(matrix3x2_assignment_operator)
{
    gil::matrix3x2<int> m1(1, 2, 3, 4, 5, 6);
    gil::matrix3x2<int> m2;
    m2 = m1;
    BOOST_TEST(m2.a == 1);
    BOOST_TEST(m2.b == 2);
    BOOST_TEST(m2.c == 3);
    BOOST_TEST(m2.d == 4);
    BOOST_TEST(m2.e == 5);
    BOOST_TEST(m2.f == 6);
}

BOOST_AUTO_TEST_CASE(matrix3x2_multiplication_assignment)
{
    gil::matrix3x2<int> m1;
    gil::matrix3x2<int> m2;
    m2 *= m1;
    BOOST_TEST(m2.a == 1);
    BOOST_TEST(m2.b == 0);
    BOOST_TEST(m2.c == 0);
    BOOST_TEST(m2.d == 1);
    BOOST_TEST(m2.e == 0);
    BOOST_TEST(m2.f == 0);

    gil::matrix3x2<int> m3(0, 0, 0, 0, 0, 0);
    m2 *= m3;
    BOOST_TEST(m2.a == 0);
    BOOST_TEST(m2.b == 0);
    BOOST_TEST(m2.c == 0);
    BOOST_TEST(m2.d == 0);
    BOOST_TEST(m2.e == 0);
    BOOST_TEST(m2.f == 0);
}

BOOST_AUTO_TEST_CASE(matrix3x2_matrix3x2_multiplication)
{
    gil::matrix3x2<int> m1;
    gil::matrix3x2<int> m2(0, 0, 0, 0, 0, 0);
    gil::matrix3x2<int> m3;
    m3 = m1 * m2;
    BOOST_TEST(m3.a == 0);
    BOOST_TEST(m3.b == 0);
    BOOST_TEST(m3.c == 0);
    BOOST_TEST(m3.d == 0);
    BOOST_TEST(m3.e == 0);
    BOOST_TEST(m3.f == 0);
}

BOOST_AUTO_TEST_CASE(matrix3x2_vector_multiplication)
{
    gil::matrix3x2<int> m1;
    gil::point<int> v1{2, 4};

    gil::point<int> v2 = v1 * m1;
    BOOST_TEST(v2.x == 2);
    BOOST_TEST(v2.y == 4);

    gil::point<int> v3 = gil::transform(m1, v1);
    BOOST_TEST(v3.x == 2);
    BOOST_TEST(v3.y == 4);
}

BOOST_AUTO_TEST_CASE(matrix3x2_get_rotate)
{
    auto m1 = gil::matrix3x2<double>::get_rotate(HALF_PI);
    BOOST_TEST(m1.a == std::cos(HALF_PI), btt::tolerance(0.03));
    BOOST_TEST(m1.b == 1);
    BOOST_TEST(m1.c == -1);
    BOOST_TEST(m1.d == std::cos(HALF_PI), btt::tolerance(0.03));
    BOOST_TEST(m1.e == 0);
    BOOST_TEST(m1.f == 0);
}

BOOST_AUTO_TEST_CASE(matrix3x2_get_scale)
{
    gil::matrix3x2<int> m1;
    m1 = gil::matrix3x2<int>::get_scale(2);
    BOOST_TEST(m1.a == 2);
    BOOST_TEST(m1.b == 0);
    BOOST_TEST(m1.c == 0);
    BOOST_TEST(m1.d == 2);
    BOOST_TEST(m1.e == 0);
    BOOST_TEST(m1.f == 0);
    m1 = gil::matrix3x2<int>::get_scale(2, 4);
    BOOST_TEST(m1.a == 2);
    BOOST_TEST(m1.d == 4);
    m1 = gil::matrix3x2<int>::get_scale(gil::point<int>{4, 8});
    BOOST_TEST(m1.a == 4);
    BOOST_TEST(m1.d == 8);
}

BOOST_AUTO_TEST_CASE(matrix3x2_get_translate)
{
    gil::matrix3x2<int> m1;
    m1 = gil::matrix3x2<int>::get_translate(2, 4);
    BOOST_TEST(m1.a == 1);
    BOOST_TEST(m1.b == 0);
    BOOST_TEST(m1.c == 0);
    BOOST_TEST(m1.d == 1);
    BOOST_TEST(m1.e == 2);
    BOOST_TEST(m1.f == 4);
    m1 = gil::matrix3x2<int>::get_translate(gil::point<int>{4, 8});
    BOOST_TEST(m1.e == 4);
    BOOST_TEST(m1.f == 8);
}

BOOST_AUTO_TEST_CASE(matrix3x2_transform)
{
    gil::matrix3x2<int> m1;
    gil::point<int> v1{2, 4};
    gil::point<int> v2 = gil::transform(m1, v1);
    BOOST_TEST(v2.x == 2);
    BOOST_TEST(v2.y == 4);
}
