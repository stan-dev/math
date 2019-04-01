#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>



TEST(MathPrimMat, test_size_type) {
  const size_type test_size_type = 1.0;
}


TEST(MathPrimMat, test_matrix_d) {
  using stan::math::matrix_d;
  matrix_d d1(3, 3);
  d1 << 1, 3, -5, 1, 3, -5, 1, 3, -5;
}


TEST(MathPrimMat, test_vector_d) {
}

TEST(MathPrimMat, test_row_vector_d) {
}
