#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>

TEST(AgradRevMatrix, sum_vector) {
  using stan::math::sum;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d(6);
  vector_v v(6);

  d << 1, 2, 3, 4, 5, 6;
  v << 1, 2, 3, 4, 5, 6;

  AVAR output;
  output = sum(d);
  EXPECT_FLOAT_EQ(21.0, output.val());

  output = sum(v);
  EXPECT_FLOAT_EQ(21.0, output.val());

  std::vector<double> grad;
  std::vector<AVAR> x(v.size());
  for (int i = 0; i < v.size(); ++i)
    x[i] = v(i);
  output.grad(x, grad);
  EXPECT_EQ(6, grad.size());
  for (int i = 0; i < 6; ++i)
    EXPECT_FLOAT_EQ(1.0, grad[i]);

  d.resize(0);
  v.resize(0);
  EXPECT_FLOAT_EQ(0.0, sum(d));
  EXPECT_FLOAT_EQ(0.0, sum(v).val());
}
TEST(AgradRevMatrix, sum_rowvector) {
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;
  using stan::math::sum;

  row_vector_d d(6);
  row_vector_v v(6);

  d << 1, 2, 3, 4, 5, 6;
  v << 1, 2, 3, 4, 5, 6;

  AVAR output;
  output = sum(d);
  EXPECT_FLOAT_EQ(21.0, output.val());

  output = sum(v);
  EXPECT_FLOAT_EQ(21.0, output.val());

  d.resize(0);
  v.resize(0);
  EXPECT_FLOAT_EQ(0.0, sum(d));
  EXPECT_FLOAT_EQ(0.0, sum(v).val());
}
TEST(AgradRevMatrix, sum_matrix) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::sum;

  matrix_d d(2, 3);
  matrix_v v(2, 3);

  d << 1, 2, 3, 4, 5, 6;
  v << 1, 2, 3, 4, 5, 6;

  AVAR output;
  output = sum(d);
  EXPECT_FLOAT_EQ(21.0, output.val());

  output = sum(v);
  EXPECT_FLOAT_EQ(21.0, output.val());

  d.resize(0, 0);
  v.resize(0, 0);
  EXPECT_FLOAT_EQ(0.0, sum(d));
  EXPECT_FLOAT_EQ(0.0, sum(v).val());
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::matrix_v m(2, 2);
  m << 1, 2, 3, 4;

  stan::math::vector_v v(3);
  v << 1, 2, 3;

  stan::math::row_vector_v rv(2);
  rv << 1, 2;

  test::check_varis_on_stack(stan::math::sum(m));
  test::check_varis_on_stack(stan::math::sum(v));
  test::check_varis_on_stack(stan::math::sum(rv));
}
