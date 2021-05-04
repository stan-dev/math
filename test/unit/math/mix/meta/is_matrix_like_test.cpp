#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <gtest/gtest.h>

TEST(MathMetaMix, is_matrix_test) {
  using stan::is_matrix;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE((is_matrix<Eigen::MatrixXd>::value));
  EXPECT_TRUE((is_matrix<Eigen::VectorXd>::value));
  EXPECT_TRUE((is_matrix<Eigen::RowVectorXd>::value));
  EXPECT_TRUE((is_matrix<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_TRUE((is_matrix<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_TRUE((is_matrix<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_TRUE((is_matrix<var_value<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_matrix<var_value<Eigen::VectorXd>>::value));
  EXPECT_TRUE((is_matrix<var_value<Eigen::RowVectorXd>>::value));
  Eigen::MatrixXd A(10, 10);
  Eigen::MatrixXd B(10, 10);
  EXPECT_TRUE((is_matrix<decltype(A * B)>::value));

  EXPECT_FALSE((is_matrix<double>::value));
  EXPECT_FALSE((is_matrix<var_value<double>>::value));
}
