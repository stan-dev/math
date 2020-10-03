#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrix, from_var_value_matrix_test) {
  Eigen::MatrixXd val(2, 3);
  val << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd adj(2, 3);
  adj << 4, 5, 6, 7, 8, 9;
  stan::math::var_value<Eigen::MatrixXd> var_value(val);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> mat_var
      = stan::math::from_var_value(var_value);
  EXPECT_MATRIX_EQ(mat_var.val(), val);
  mat_var.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(var_value.adj(), adj);
}

TEST(AgradRevMatrix, from_var_value_vector_test) {
  Eigen::VectorXd val(3);
  val << 1, 2, 3;
  Eigen::VectorXd adj(3);
  adj << 4, 5, 6;
  stan::math::var_value<Eigen::VectorXd> var_value(val);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> mat_var
      = stan::math::from_var_value(var_value);
  EXPECT_MATRIX_EQ(mat_var.val(), val);
  mat_var.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(var_value.adj(), adj);
}

TEST(AgradRevMatrix, from_var_value_row_vector_test) {
  Eigen::RowVectorXd val(3);
  val << 1, 2, 3;
  Eigen::RowVectorXd adj(3);
  adj << 7, 8, 9;
  stan::math::var_value<Eigen::RowVectorXd> var_value(val);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> mat_var
      = stan::math::from_var_value(var_value);
  EXPECT_MATRIX_EQ(mat_var.val(), val);
  mat_var.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(var_value.adj(), adj);
}
