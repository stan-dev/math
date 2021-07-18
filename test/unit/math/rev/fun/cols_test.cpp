#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(AgradRevMatrix, cols_vector) {
  using stan::math::cols;
  using stan::math::row_vector_v;
  using stan::math::var_value;
  using stan::math::vector_v;

  Eigen::VectorXd vec(5);
  vec << 0, 1, 2, 3, 4;
  vector_v vec_v(vec);
  var_value<Eigen::VectorXd> vec_vv(vec);
  EXPECT_EQ(1U, cols(vec_v));
  EXPECT_EQ(1U, cols(vec_vv));

  vec_v.resize(0);
  EXPECT_EQ(1U, cols(vec_v));
}
TEST(AgradRevMatrix, cols_rowvector) {
  using stan::math::cols;
  using stan::math::row_vector_v;
  using stan::math::var_value;
  Eigen::RowVectorXd r(5);
  r << 0, 1, 2, 3, 4;
  row_vector_v rv(r);
  var_value<Eigen::RowVectorXd> r_vv(r);
  EXPECT_EQ(5U, cols(rv));
  EXPECT_EQ(5U, cols(r_vv));

  rv.resize(0);
  EXPECT_EQ(0U, cols(rv));
}

TEST(AgradRevMatrix, cols_matrix) {
  using stan::math::cols;
  using stan::math::matrix_v;
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;
  matrix_v m_v(m);
  EXPECT_EQ(3U, cols(m_v));
  var_value<Eigen::MatrixXd> m_vv(m);
  EXPECT_EQ(3U, cols(m_vv));

  m.resize(5, 0);
  EXPECT_EQ(0U, cols(m));
}
