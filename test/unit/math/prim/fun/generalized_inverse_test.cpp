#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, Zero) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d ginv = stan::math::generalized_inverse(m0);
  EXPECT_EQ(0, ginv.rows());
  EXPECT_EQ(0, ginv.cols());
}

TEST(MathMatrixPrim, Equal) {
  using stan::math::generalized_inverse;

  stan::math::matrix_d m1(3, 2);
  m1 << 1, 2, 2, 4, 1, 2;

  stan::math::matrix_d m2(2, 3);
  m2 << 1/30, 1/15, 1/30, 1/15, 2/15, 1/15;
  EXPECT_MATRIX_FLOAT_EQ(generalized_inverse(m1), m2);
}
