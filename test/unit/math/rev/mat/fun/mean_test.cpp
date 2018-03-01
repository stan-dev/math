#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, mean_vector) {
  using stan::math::mean;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d1(3);
  vector_v v1(3);

  d1 << 100, 0, -3;
  v1 << 100, 0, -3;

  AVAR output;
  output = mean(d1);
  EXPECT_FLOAT_EQ(97.0 / 3.0, output.val());

  output = mean(v1);
  EXPECT_FLOAT_EQ(97.0 / 3.0, output.val());
}
TEST(AgradRevMatrix, mean_vector_exception) {
  using stan::math::mean;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d;
  vector_v v;
  EXPECT_THROW(mean(d), std::invalid_argument);
  EXPECT_THROW(mean(v), std::invalid_argument);
}
TEST(AgradRevMatrix, mean_rowvector) {
  using stan::math::mean;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d d1(3);
  row_vector_v v1(3);

  d1 << 100, 0, -3;
  v1 << 100, 0, -3;

  AVAR output;
  output = mean(d1);
  EXPECT_FLOAT_EQ(97.0 / 3.0, output.val());

  output = mean(v1);
  EXPECT_FLOAT_EQ(97.0 / 3.0, output.val());
}
TEST(AgradRevMatrix, mean_rowvector_exception) {
  using stan::math::mean;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d d;
  row_vector_v v;
  EXPECT_THROW(mean(d), std::invalid_argument);
  EXPECT_THROW(mean(v), std::invalid_argument);
}
TEST(AgradRevMatrix, mean_matrix) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mean;

  matrix_d d1(3, 1);
  matrix_v v1(1, 3);

  d1 << 100, 0, -3;
  v1 << 100, 0, -3;

  AVAR output;
  output = mean(d1);
  EXPECT_FLOAT_EQ(97.0 / 3.0, output.val());

  output = mean(v1);
  EXPECT_FLOAT_EQ(97.0 / 3.0, output.val());
}
TEST(AgradRevMatrix, mean_matrix_exception) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mean;

  matrix_d d;
  matrix_v v;
  EXPECT_THROW(mean(d), std::invalid_argument);
  EXPECT_THROW(mean(v), std::invalid_argument);
}
TEST(AgradRevMatrix, meanStdVector) {
  // should use arg-dep lookup
  using stan::math::mean;
  AVEC x(0);
  EXPECT_THROW(mean(x), std::invalid_argument);
  x.push_back(1.0);
  EXPECT_FLOAT_EQ(1.0, mean(x).val());
  x.push_back(2.0);
  EXPECT_FLOAT_EQ(1.5, mean(x).val());

  AVEC y = createAVEC(1.0, 2.0);
  AVAR f = mean(y);
  VEC grad = cgrad(f, y[0], y[1]);
  EXPECT_FLOAT_EQ(0.5, grad[0]);
  EXPECT_FLOAT_EQ(0.5, grad[1]);
  EXPECT_EQ(2U, grad.size());
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::vector_v v(3);
  v << -100, 0, 1;
  stan::math::row_vector_v rv(3);
  rv << -100, 0, 1;
  stan::math::matrix_v m(2, 3);
  m << -100, 0, 1, 20, -40, 2;
  test::check_varis_on_stack(stan::math::mean(v));
  test::check_varis_on_stack(stan::math::mean(rv));
  test::check_varis_on_stack(stan::math::mean(m));
}
