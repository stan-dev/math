#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, inverse_exception) {
  using stan::math::inverse;

  stan::math::matrix_d m0(0, 0);
  EXPECT_THROW(inverse(m0), std::invalid_argument);

  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(inverse(m1), std::invalid_argument);
}
