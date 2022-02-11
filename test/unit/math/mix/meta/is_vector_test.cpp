#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <gtest/gtest.h>

TEST(MathMetaMix, is_row_vector_test) {
  using stan::is_row_vector;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE((is_row_vector<Eigen::RowVectorXd>::value));
  EXPECT_TRUE((is_row_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_TRUE((is_row_vector<var_value<Eigen::RowVectorXd>>::value));
  Eigen::MatrixXd A(10, 10);
  Eigen::MatrixXd B(10, 10);
  EXPECT_FALSE((is_row_vector<decltype(A * B)>::value));
  EXPECT_FALSE((is_row_vector<var_value<Eigen::MatrixXd>>::value));
  EXPECT_FALSE((is_row_vector<var_value<Eigen::VectorXd>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::MatrixXd>::value));
  EXPECT_FALSE((is_row_vector<Eigen::VectorXd>::value));

  EXPECT_FALSE((is_row_vector<double>::value));
  EXPECT_FALSE((is_row_vector<var_value<double>>::value));
}

TEST(MathMetaMix, is_col_vector_test) {
  using stan::is_col_vector;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE((is_col_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_TRUE((is_col_vector<Eigen::VectorXd>::value));
  EXPECT_TRUE((is_col_vector<var_value<Eigen::VectorXd>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::RowVectorXd>::value));
  EXPECT_FALSE((is_col_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_FALSE((is_col_vector<var_value<Eigen::RowVectorXd>>::value));
  Eigen::MatrixXd A(10, 10);
  Eigen::MatrixXd B(10, 10);
  EXPECT_FALSE((is_col_vector<decltype(A * B)>::value));
  EXPECT_FALSE((is_col_vector<var_value<Eigen::MatrixXd>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::MatrixXd>::value));

  EXPECT_FALSE((is_col_vector<double>::value));
  EXPECT_FALSE((is_col_vector<var_value<double>>::value));
}
