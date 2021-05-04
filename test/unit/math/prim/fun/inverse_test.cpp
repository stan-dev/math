#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, inverse) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d inv = stan::math::inverse(m0);
  EXPECT_EQ(0, inv.rows());
  EXPECT_EQ(0, inv.cols());
}

TEST(MathMatrixPrim, inverse_exception) {
  using stan::math::inverse;

  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(inverse(m1), std::invalid_argument);
}
