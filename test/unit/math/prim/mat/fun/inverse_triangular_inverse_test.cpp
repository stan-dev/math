#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, inverse_exception) {
  stan::math::matrix_d m1(2,3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_gpu m2(m1);
  using stan::math::lower_triangular_inverse;
  EXPECT_THROW(lower_triangular_inverse(m2),std::invalid_argument);
}
