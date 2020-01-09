#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, diag_matrix_exception) {
  stan::math::vector_d v0;

  using stan::math::diag_matrix;
  EXPECT_NO_THROW(diag_matrix(v0));
}
