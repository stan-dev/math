#include <stan/math/mix.hpp>
#include <test/unit/math/mix/util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST_F(mathMix, ffd_vector) {
  using stan::math::row_vector_ffd;
  using stan::math::rows;
  using stan::math::vector_ffd;

  vector_ffd v(5);
  v << 0, 1, 2, 3, 4;
  EXPECT_EQ(5U, rows(v));

  v.resize(0);
  EXPECT_EQ(0U, rows(v));
}
TEST_F(mathMix, ffd_rowvector) {
  using stan::math::row_vector_ffd;
  using stan::math::rows;

  row_vector_ffd rv(5);
  rv << 0, 1, 2, 3, 4;
  EXPECT_EQ(1U, rows(rv));

  rv.resize(0);
  EXPECT_EQ(1U, rows(rv));
}
TEST_F(mathMix, ffd_matrix) {
  using stan::math::matrix_ffd;
  using stan::math::rows;

  matrix_ffd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;
  EXPECT_EQ(2U, rows(m));

  m.resize(0, 2);
  EXPECT_EQ(0U, rows(m));
}
