#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrix, col_varmat) {
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::MatrixXd> m_v(m);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::col(m_v, 2).val(), stan::math::col(m, 2));
}

TEST(AgradRevMatrix, row_varmat) {
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::MatrixXd> m_v(m);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::row(m_v, 2).val(), stan::math::row(m, 2));
}

TEST(AgradRevMatrix, dims_varmat) {
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::MatrixXd> m_v(m);
  EXPECT_STD_VECTOR_EQ(stan::math::dims(m_v), stan::math::dims(m));
}
TEST(AgradRevMatrix, num_elements_varmat) {
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::MatrixXd> m_v(m);
  EXPECT_EQ(stan::math::num_elements(m_v), stan::math::num_elements(m));
}

TEST(AgradRevMatrix, diagonal_varmat) {
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::MatrixXd> m_v(m);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::diagonal(m_v).val(),
                         stan::math::diagonal(m));
}

TEST(AgradRevMatrix, subsets_varmat) {
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::MatrixXd> m_v(m);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::sub_row(m_v, 1, 1, 2).val(),
                         stan::math::sub_row(m, 1, 1, 2));
  EXPECT_MATRIX_FLOAT_EQ(stan::math::sub_col(m_v, 1, 1, 2).val(),
                         stan::math::sub_col(m, 1, 1, 2));
}

TEST(AgradRevMatrix, transpose_varmat) {
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::MatrixXd> m_v(m);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::transpose(m_v).val(),
                         stan::math::transpose(m));
}

TEST(AgradRevMatrix, head_varmat) {
  using stan::math::var_value;
  Eigen::VectorXd vec(6);
  vec << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::VectorXd> vec_v(vec);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::head(vec_v, 2).val(),
                         stan::math::head(vec, 2));
}

TEST(AgradRevMatrix, tail_varmat) {
  using stan::math::var_value;
  Eigen::VectorXd vec(6);
  vec << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::VectorXd> vec_v(vec);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::tail(vec_v, 2).val(),
                         stan::math::tail(vec, 2));
}

TEST(AgradRevMatrix, reverse_varmat) {
  using stan::math::var_value;
  Eigen::VectorXd vec(6);
  vec << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::VectorXd> vec_v(vec);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::reverse(vec_v).val(),
                         stan::math::reverse(vec));
}

TEST(AgradRevMatrix, segment_varmat) {
  using stan::math::var_value;
  Eigen::VectorXd vec(6);
  vec << 0, 1, 2, 3, 4, 5;

  var_value<Eigen::VectorXd> vec_v(vec);
  EXPECT_MATRIX_FLOAT_EQ(stan::math::segment(vec_v, 2, 2).val(),
                         stan::math::segment(vec, 2, 2));
}
