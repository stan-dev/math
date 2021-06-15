#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(AgradRevMatrix, rows_vector) {
  using stan::math::row_vector_v;
  using stan::math::rows;
  using stan::math::var_value;
  using stan::math::vector_v;
  Eigen::VectorXd vec(5);
  vec << 0, 1, 2, 3, 4;
  vector_v vec_v(5);
  var_value<Eigen::VectorXd> vec_vv(vec);
  EXPECT_EQ(5U, rows(vec_v));
  EXPECT_EQ(5U, rows(vec_vv));

  vec_v.resize(0);
  EXPECT_EQ(0U, rows(vec_v));
}
TEST(AgradRevMatrix, rows_rowvector) {
  using stan::math::row_vector_v;
  using stan::math::rows;
  using stan::math::var_value;

  Eigen::RowVectorXd r(5);
  r << 0, 1, 2, 3, 4;
  row_vector_v r_v(r);
  var_value<Eigen::RowVectorXd> r_vv(r);
  EXPECT_EQ(1U, rows(r_v));
  EXPECT_EQ(1U, rows(r_v));

  r_v.resize(0);
  EXPECT_EQ(1U, rows(r_v));
}

TEST(AgradRevMatrix, rows_matrix) {
  using stan::math::matrix_v;
  using stan::math::rows;
  using stan::math::var_value;
  Eigen::MatrixXd m(2, 3);
  m << 0, 1, 2, 3, 4, 5;
  matrix_v m_v(2, 3);
  EXPECT_EQ(2U, rows(m_v));

  m_v.resize(0, 2);
  EXPECT_EQ(0U, rows(m_v));
}
