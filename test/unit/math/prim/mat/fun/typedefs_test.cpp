#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>

TEST(MathPrimMat, test_size_type) {
  const size_type test_size_type = 1.0;
}

TEST(MathPrimMat, test_matrix_d) {
  using stan::math::matrix_d;
  matrix_d d(3, 3);
  d << 1, 3, -5, 1, 3, -5, 1, 3, -5;
}

TEST(MathPrimMat, test_vector_d) {
  using stan::math::vector_d;
  vector_d d(3);
  d << 1, 2, 3;
}

TEST(MathPrimMat, test_row_vector_d) {
  using stan::math::row_vector_d;
  row_vector_d d(4);
  d << 1, 2, 3, 5;
}
