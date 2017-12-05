#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, transpose) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;

  using stan::math::transpose;
  EXPECT_NO_THROW(transpose(v0));
  EXPECT_NO_THROW(transpose(rv0));
  EXPECT_NO_THROW(transpose(m0));
}
