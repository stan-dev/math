#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMatrixPrimMat, svd_V) {
    using stan::math::svd_V;
    using stan::math::matrix_d;
    matrix_d m23(2, 3);
    m23 << 1, 3, -5, 7, 9, -11;

    matrix_d m23_V(3, 3);
    m23_V << -0.4117624, 0.8147301, -0.4082483, -0.5638395,
        0.1241705, 0.8164966, 0.7159167, 0.5663891, 0.4082483; //Generated using R base::svd
    EXPECT_MATRIX_FLOAT_EQ(m23_V, svd_V(m23));
}