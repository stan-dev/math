#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, inverse_ldlt_val) {
  stan::math::matrix_d Ad(2, 2);
  stan::math::matrix_d I;

  Ad << 2.0, 3.0, 3.0, 7.0;

  auto ldlt_Ad = stan::math::make_ldlt_factor(Ad);

  I = inverse_ldlt(ldlt_Ad);
  EXPECT_NEAR(1.4, I(0, 0), 1.0E-12);
  EXPECT_NEAR(-0.6, I(0, 1), 1.0E-12);
  EXPECT_NEAR(-0.6, I(1, 0), 1.0E-12);
  EXPECT_NEAR(0.4, I(1, 1), 1.0e-12);
}

TEST(MathMatrixPrimMat, inverse_ldlt_val_0x0) {
  stan::math::matrix_d A(0, 0);

  auto ldlt_A = stan::math::make_ldlt_factor(A);

  auto M = inverse_ldlt(ldlt_A);
  EXPECT_EQ(0, M.rows());
  EXPECT_EQ(0, M.cols());
}
