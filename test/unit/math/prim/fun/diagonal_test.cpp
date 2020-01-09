#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, diagonal) {
  stan::math::matrix_d m0;

  using stan::math::diagonal;
  EXPECT_NO_THROW(diagonal(m0));
}
