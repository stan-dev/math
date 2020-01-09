#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, singular_values) {
  using stan::math::singular_values;

  stan::math::matrix_d m0(0, 0);
  EXPECT_NO_THROW(singular_values(m0));

  stan::math::matrix_d m1(1, 1);
  m1 << 1.0;
  EXPECT_NO_THROW(singular_values(m1));
}
