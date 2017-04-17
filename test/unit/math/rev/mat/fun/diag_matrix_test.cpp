#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix,diagMatrix) {
  using stan::math::diag_matrix;
  using stan::math::matrix_v;
  using stan::math::vector_d;
  using stan::math::vector_v;

  EXPECT_EQ(0,diag_matrix(vector_v()).size());
  EXPECT_EQ(4,diag_matrix(vector_v(2)).size());
  EXPECT_EQ(0,diag_matrix(vector_d()).size());
  EXPECT_EQ(4,diag_matrix(vector_d(2)).size());

  vector_v v(3);
  v << 1, 4, 9;
  matrix_v m = diag_matrix(v);
  EXPECT_EQ(1,m(0,0).val());
  EXPECT_EQ(4,m(1,1).val());
  EXPECT_EQ(9,m(2,2).val());
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::vector_v v(3);
  v << 1, 4, 9;
  test::check_varis_on_stack(stan::math::diag_matrix(v));
}
