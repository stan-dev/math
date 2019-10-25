#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, inverse_exception) {
  stan::math::vector_d v0;

  using stan::math::diag_matrix;
  EXPECT_NO_THROW(diag_matrix(v0));
}
