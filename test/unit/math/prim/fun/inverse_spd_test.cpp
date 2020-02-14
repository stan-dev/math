#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, inverse_spd) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d inv = stan::math::inverse_spd(m0);
  EXPECT_EQ(0, inv.rows());
  EXPECT_EQ(0, inv.cols());
}

TEST(MathMatrixPrim, inverse_spd_exception) {
  using stan::math::inverse_spd;

  stan::math::matrix_d m1(2, 3);

  // non-square
  m1 << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(inverse_spd(m1), std::invalid_argument);

  stan::math::matrix_d m2(3, 3);

  // non-symmetric
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  EXPECT_THROW(inverse_spd(m2), std::domain_error);

  // not positive definite
  m2 << 1, 2, 3, 2, 4, 5, 3, 5, 6;
  EXPECT_THROW(inverse_spd(m2), std::domain_error);
}
