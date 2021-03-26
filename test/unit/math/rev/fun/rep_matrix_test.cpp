#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <limits>
#include <string>
#include <vector>

TEST(MathMixMatFun, repMatrixVar) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  auto x_var = var(1.0);
  Eigen::Matrix<var, Eigen::Dynamic, 1> x_vec_var = Eigen::VectorXd(2);
  Eigen::Matrix<var, 1, Eigen::Dynamic> x_row_vec_var = Eigen::RowVectorXd(2);
  auto x1 = rep_matrix(x_var, 2, 3);
  auto x2 = rep_matrix<Eigen::Matrix<var, -1, -1>>(x_var, 2, 3);
  auto x3 = rep_matrix(x_vec_var, 3);
  auto x4 = rep_matrix(x_row_vec_var, 3);
  EXPECT_TRUE(stan::is_eigen<decltype(x1)>::value);
  EXPECT_TRUE(stan::is_eigen<decltype(x2)>::value);
  EXPECT_TRUE(stan::is_eigen<decltype(x3)>::value);
  EXPECT_TRUE(stan::is_eigen<decltype(x4)>::value);

  EXPECT_EQ(x1.rows(), 2);
  EXPECT_EQ(x1.cols(), 3);
  EXPECT_EQ(x2.rows(), 2);
  EXPECT_EQ(x2.cols(), 3);
  EXPECT_EQ(x3.rows(), 2);
  EXPECT_EQ(x3.cols(), 3);
  EXPECT_EQ(x4.rows(), 3);
  EXPECT_EQ(x4.cols(), 2);
}

TEST(MathMixMatFun, repVarMatrix) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  auto x_var = var(1.0);
  auto x = rep_matrix<var_value<Eigen::MatrixXd>>(x_var, 5, 5);
  EXPECT_TRUE(stan::is_var_matrix<decltype(x)>::value);
  auto x_sum = sum(x);
  x_sum.grad();

  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_EQ(x_var.adj(), 25.0);
}

TEST(MathMixMatFun, repVarMatrixVec) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  var_value<Eigen::VectorXd> x_var(Eigen::VectorXd::Ones(5));
  auto x = rep_matrix(x_var, 5);
  EXPECT_TRUE(stan::is_var_matrix<decltype(x)>::value);
  auto x_sum = sum(x);
  x_sum.grad();
  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_MATRIX_EQ(x_var.val(), Eigen::VectorXd::Ones(5));
  Eigen::VectorXd expected_x_var_adjs(5);
  expected_x_var_adjs << 5, 5, 5, 5, 5;
  EXPECT_MATRIX_EQ(x_var.adj(), expected_x_var_adjs);
}

TEST(MathMixMatFun, repVarMatrixRowVec) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  var_value<Eigen::RowVectorXd> x_var(Eigen::RowVectorXd::Ones(5));
  auto x = rep_matrix(x_var, 5);
  EXPECT_TRUE(stan::is_var_matrix<decltype(x)>::value);
  auto x_sum = sum(x);
  x_sum.grad();

  EXPECT_EQ(x_sum.val(), 25.0);
  Eigen::RowVectorXd expected_x_var_adjs(5);
  expected_x_var_adjs << 5, 5, 5, 5, 5;
  EXPECT_MATRIX_EQ(x_var.adj(), expected_x_var_adjs);
}
