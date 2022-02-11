#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, plus) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;

  EXPECT_EQ(0, stan::math::plus(v0).size());
  EXPECT_EQ(0, stan::math::plus(rv0).size());
  EXPECT_EQ(0, stan::math::plus(m0).size());
}
