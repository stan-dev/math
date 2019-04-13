#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, mdivide_right_tri_low_val1) {
  using stan::math::mdivide_right_tri_low;
  stan::math::matrix_d Ad(2, 2);
  stan::math::matrix_d I;

  Ad << 2.0, 0.0, 5.0, 7.0;

  I = mdivide_right_tri_low(Ad, Ad);
  EXPECT_EQ(1, I(0, 0));
  EXPECT_EQ(0, I(0, 1));
  EXPECT_EQ(0, I(1, 0));
  EXPECT_EQ(1, I(1, 1));
}
