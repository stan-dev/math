#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, elt_divide_vec_vv) {
  using stan::math::elt_divide;
  using stan::math::vector_v;

  vector_v x(2);
  x << 2, 5;
  vector_v y(2);
  y << 10, 100;
  AVEC x_ind = createAVEC(x(0), x(1), y(0), y(1));
  vector_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0).val());
  EXPECT_FLOAT_EQ(0.05, z(1).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(1.0/10.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(2.0 / (- 10.0 * 10.0), g[2]);
  EXPECT_FLOAT_EQ(0.0, g[3]);
}
TEST(AgradRevMatrix, elt_divide_vec_vd) {
  using stan::math::elt_divide;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_v x(2);
  x << 2, 5;
  vector_d y(2);
  y << 10, 100;
  AVEC x_ind = createAVEC(x(0), x(1));
  vector_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0).val());
  EXPECT_FLOAT_EQ(0.05, z(1).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(1.0/10.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
}
TEST(AgradRevMatrix, elt_divide_vec_dv) {
  using stan::math::elt_divide;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d x(2);
  x << 2, 5;
  vector_v y(2);
  y << 10, 100;
  AVEC x_ind = createAVEC(y(0), y(1));
  vector_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0).val());
  EXPECT_FLOAT_EQ(0.05, z(1).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(2.0 / (- 10.0 * 10.0), g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
}

TEST(AgradRevMatrix, elt_divide_rowvec_vv) {
  using stan::math::elt_divide;
  using stan::math::row_vector_v;

  row_vector_v x(2);
  x << 2, 5;
  row_vector_v y(2);
  y << 10, 100;
  AVEC x_ind = createAVEC(x(0), x(1), y(0), y(1));
  row_vector_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0).val());
  EXPECT_FLOAT_EQ(0.05, z(1).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(1.0/10.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(2.0 / (- 10.0 * 10.0), g[2]);
  EXPECT_FLOAT_EQ(0.0, g[3]);
}
TEST(AgradRevMatrix, elt_divide_rowvec_vd) {
  using stan::math::elt_divide;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_v x(2);
  x << 2, 5;
  row_vector_d y(2);
  y << 10, 100;
  AVEC x_ind = createAVEC(x(0), x(1));
  row_vector_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0).val());
  EXPECT_FLOAT_EQ(0.05, z(1).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(1.0/10.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
}
TEST(AgradRevMatrix, elt_divide_rowvec_dv) {
  using stan::math::elt_divide;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d x(2);
  x << 2, 5;
  row_vector_v y(2);
  y << 10, 100;
  AVEC x_ind = createAVEC(y(0), y(1));
  row_vector_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0).val());
  EXPECT_FLOAT_EQ(0.05, z(1).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(2.0 / (- 10.0 * 10.0), g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
}


TEST(AgradRevMatrix, elt_divide_mat_vv) {
  using stan::math::elt_divide;
  using stan::math::matrix_v;

  matrix_v x(2, 3);
  x << 2, 5, 7, 13, 29, 112;
  matrix_v y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;
  AVEC x_ind = createAVEC(x(0, 0), x(0, 1), y(0, 0), y(0, 1));
  matrix_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0, 0).val());
  EXPECT_FLOAT_EQ(0.05, z(0, 1).val());
  EXPECT_FLOAT_EQ(112.0/1000000.0, z(1, 2).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(1.0/10.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(2.0 / (- 10.0 * 10.0), g[2]);
  EXPECT_FLOAT_EQ(0.0, g[3]);
}
TEST(AgradRevMatrix, elt_divide_mat_vd) {
  using stan::math::elt_divide;
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_v x(2, 3);
  x << 2, 5, 7, 13, 29, 112;
  matrix_d y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;
  AVEC x_ind = createAVEC(x(0, 0), x(0, 1));
  matrix_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0, 0).val());
  EXPECT_FLOAT_EQ(0.05, z(0, 1).val());
  EXPECT_FLOAT_EQ(112.0/1000000.0, z(1, 2).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(1.0/10.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
}
TEST(AgradRevMatrix, elt_divide_mat_dv) {
  using stan::math::elt_divide;
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_d x(2, 3);
  x << 2, 5, 7, 13, 29, 112;
  matrix_v y(2, 3);
  y << 10, 100, 1000, 10000, 100000, 1000000;
  AVEC x_ind = createAVEC(y(0, 0), y(0, 1));
  matrix_v z = elt_divide(x, y);
  EXPECT_FLOAT_EQ(0.2, z(0, 0).val());
  EXPECT_FLOAT_EQ(0.05, z(0, 1).val());
  EXPECT_FLOAT_EQ(112.0/1000000.0, z(1, 2).val());

  VEC g = cgradvec(z(0), x_ind);
  EXPECT_FLOAT_EQ(2.0 / (- 10.0 * 10.0), g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
}
TEST(AgradRevMatrix, elt_divide_mat_scal_dv) {
  using stan::math::elt_divide;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::var;

  matrix_d x(2, 3);
  x << 2, 5, 7, 13, 29, 112;

  var y = 10;

  matrix_v z = elt_divide(x, y);
  z.sum().grad();
  EXPECT_FLOAT_EQ(x.sum() * (-1.0 / 100),  y.adj());
}
TEST(AgradRevMatrix, elt_divide_vec_scal_dv) {
  using stan::math::elt_divide;
  using stan::math::vector_d;
  using stan::math::vector_v;
  using stan::math::var;

  vector_d x(6);
  x << 2, 5, 7, 13, 29, 112;

  var y = 10;

  vector_v z = elt_divide(x, y);
  z.sum().grad();
  EXPECT_FLOAT_EQ(x.sum() * (-1.0 / 100),  y.adj());
}
TEST(AgradRevMatrix, elt_divide_row_vec_scal_dv) {
  using stan::math::elt_divide;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;
  using stan::math::var;

  row_vector_d x(6);
  x << 2, 5, 7, 13, 29, 112;

  var y = 10;

  row_vector_v z = elt_divide(x, y);
  z.sum().grad();
  EXPECT_FLOAT_EQ(x.sum() * (-1.0 / 100),  y.adj());
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::value_of;
  using stan::math::elt_divide;
  using stan::math::vector_v;
  vector_v x(2);
  x << 2, 5;
  vector_v y(2);
  y << 10, 100;

  test::check_varis_on_stack(stan::math::elt_divide(x, y));
  test::check_varis_on_stack(stan::math::elt_divide(x, value_of(y)));
  test::check_varis_on_stack(stan::math::elt_divide(value_of(x), y));
}
