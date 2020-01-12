#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, inverse_spd_exception) {
  using stan::math::inverse_spd;

  // empty
  stan::math::matrix_d m0(0, 0);
  EXPECT_THROW(inverse_spd(m0), std::invalid_argument);

  stan::math::matrix_d m1(2, 3);

  // non-square
  m1 << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(inverse_spd(m1), std::invalid_argument);

  stan::math::matrix_d m2(3, 3);

  // non-symmetric
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  EXPECT_THROW(inverse_spd(m1), std::invalid_argument);

  // not positive definite
  m2 << 1, 2, 3, 2, 4, 5, 3, 5, 6;
  EXPECT_THROW(inverse_spd(m1), std::invalid_argument);
}
