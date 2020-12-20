#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMatrixPrimMat, svd_U) {
    using stan::math::svd_U;
    using stan::math::matrix_d;
    matrix_d m23(2, 3);
    m23 << 1, 3, -5, 7, 9, -11;

    matrix_d m23_U(2, 2);
    m23_U << -0.3378432, -0.9412024, -0.9412024, 0.3378432; //Generated using R base::svd
    EXPECT_MATRIX_FLOAT_EQ(m23_U, svd_U(m23));
}