#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>

TEST(AgradRevMatrix, cols_vector) {
  using stan::math::vector_v;
  using stan::math::row_vector_v;
  using stan::math::cols;

  vector_v v(5);
  v << 0, 1, 2, 3, 4;
  EXPECT_EQ(1U, cols(v));

  v.resize(0);
  EXPECT_EQ(1U, cols(v));
}
TEST(AgradRevMatrix, cols_rowvector) {
  using stan::math::row_vector_v;
  using stan::math::cols;

  row_vector_v rv(5);
  rv << 0, 1, 2, 3, 4;
  EXPECT_EQ(5U, cols(rv));

  rv.resize(0);
  EXPECT_EQ(0U, cols(rv));
}
TEST(AgradRevMatrix, cols_matrix) {
  using stan::math::matrix_v;
  using stan::math::cols;

  matrix_v m(2, 3);
  m << 0, 1, 2, 3, 4, 5;
  EXPECT_EQ(3U, cols(m));

  m.resize(5, 0);
  EXPECT_EQ(0U, cols(m));
}
