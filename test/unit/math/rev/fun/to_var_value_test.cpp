#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrix, to_var_value_matrix_test) {
  Eigen::MatrixXd val(2, 3);
  val << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd adj(2, 3);
  val << 4, 5, 6, 7, 8, 9;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> mat_var = val;
  stan::math::var_value<Eigen::MatrixXd> var_value
      = stan::math::to_var_value(mat_var);
  EXPECT_MATRIX_EQ(var_value.val(), val);
  var_value.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(mat_var.adj(), adj);
}

TEST(AgradRevMatrix, to_var_value_vector_test) {
  Eigen::VectorXd val(3);
  val << 1, 2, 3;
  Eigen::VectorXd adj(3);
  val << 7, 8, 9;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> mat_var = val;
  stan::math::var_value<Eigen::VectorXd> var_value
      = stan::math::to_var_value(mat_var);
  EXPECT_MATRIX_EQ(var_value.val(), val);
  var_value.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(mat_var.adj(), adj);
}

TEST(AgradRevMatrix, to_var_value_row_vector_test) {
  Eigen::RowVectorXd val(3);
  val << 1, 2, 3;
  Eigen::RowVectorXd adj(3);
  val << 7, 8, 9;
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> mat_var = val;
  stan::math::var_value<Eigen::RowVectorXd> var_value
      = stan::math::to_var_value(mat_var);
  EXPECT_MATRIX_EQ(var_value.val(), val);
  var_value.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(mat_var.adj(), adj);
}
