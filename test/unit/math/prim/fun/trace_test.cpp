#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, trace) {
  using stan::math::trace;
  stan::math::matrix_d m;
  EXPECT_FLOAT_EQ(0.0, trace(m));
  m = stan::math::matrix_d(1, 1);
  m << 2.3;
  EXPECT_FLOAT_EQ(2.3, trace(m));
  m = stan::math::matrix_d(2, 3);
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_FLOAT_EQ(6.0, trace(m));
}
