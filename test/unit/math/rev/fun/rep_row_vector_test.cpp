#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMixMatFun, repRowVectorVar) {
  using stan::math::rep_row_vector;
  using stan::math::sum;
  using stan::math::var;
  var x_var = var(1.0);
  Eigen::Matrix<var, 1, -1> x1 = rep_row_vector(x_var, 3);
  Eigen::Matrix<var, 1, -1> x2
      = rep_row_vector<Eigen::Matrix<var, 1, -1>>(x_var, 3);
  EXPECT_TRUE(stan::is_eigen<decltype(x1)>::value);
  EXPECT_TRUE(stan::is_eigen<decltype(x2)>::value);

  EXPECT_EQ(x1.rows(), 1);
  EXPECT_EQ(x2.rows(), 1);
  EXPECT_EQ(x1.cols(), 3);
  EXPECT_EQ(x2.cols(), 3);
}

TEST(MathMixMatFun, repRowVarVector) {
  using stan::math::rep_row_vector;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  var x_var = var(1.0);
  var_value<Eigen::RowVectorXd> x
      = rep_row_vector<var_value<Eigen::RowVectorXd>>(x_var, 5);
  EXPECT_TRUE(stan::is_var_matrix<decltype(x)>::value);
  auto x_sum = sum(x);
  x_sum.grad();

  EXPECT_EQ(x_sum.val(), 5.0);
  EXPECT_EQ(x_sum.adj(), 1.0);
  EXPECT_EQ(x_var.val(), 1.0);
  EXPECT_EQ(x_var.adj(), 5.0);
}
