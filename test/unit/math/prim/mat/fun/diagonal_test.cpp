#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>

TEST(MathMatrix, diagonal) {
  stan::math::matrix_d m0;

  using stan::math::diagonal;
  EXPECT_NO_THROW(diagonal(m0));
}
