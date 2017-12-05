#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, col_v) {
  using stan::math::col;
  using stan::math::matrix_v;
  using stan::math::vector_v;

  matrix_v y(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  vector_v z = col(y, 1);
  EXPECT_EQ(2, z.size());
  EXPECT_FLOAT_EQ(1.0, z[0].val());
  EXPECT_FLOAT_EQ(4.0, z[1].val());

  vector_v w = col(y, 2);
  EXPECT_EQ(2, w.size());
  EXPECT_EQ(2.0, w[0].val());
  EXPECT_EQ(5.0, w[1].val());
}
TEST(AgradRevMatrix, col_v_exc0) {
  using stan::math::col;
  using stan::math::matrix_v;

  matrix_v y(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(col(y, 0), std::out_of_range);
  EXPECT_THROW(col(y, 7), std::out_of_range);
}
TEST(AgradRevMatrix, col_v_excHigh) {
  using stan::math::col;
  using stan::math::matrix_v;

  matrix_v y(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(col(y, 0), std::out_of_range);
  EXPECT_THROW(col(y, 5), std::out_of_range);
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::matrix_v y(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  test::check_varis_on_stack(stan::math::col(y, 2));
}
