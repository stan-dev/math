#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, mdivide_left_ldlt_val) {
  stan::math::matrix_d Ad(2, 2);
  stan::math::matrix_d I;

  Ad << 2.0, 3.0, 3.0, 7.0;

  auto ldlt_Ad = stan::math::make_ldlt_factor(Ad);

  I = mdivide_left_ldlt(ldlt_Ad, Ad);
  EXPECT_NEAR(1.0, I(0, 0), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1), 1.0e-12);
}

TEST(MathMatrixPrimMat, mdivide_left_ldlt_val_0x0) {
  stan::math::matrix_d A(0, 0), B(0, 0), C(0, 2);

  auto ldlt_A = stan::math::make_ldlt_factor(A);

  auto M = mdivide_left_ldlt(ldlt_A, B);
  EXPECT_EQ(0, M.rows());
  EXPECT_EQ(B.cols(), M.cols());

  auto N = mdivide_left_ldlt(ldlt_A, C);
  EXPECT_EQ(0, N.rows());
  EXPECT_EQ(C.cols(), N.cols());
}
