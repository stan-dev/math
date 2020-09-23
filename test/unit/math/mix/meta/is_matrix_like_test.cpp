#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <gtest/gtest.h>

TEST(MathMetaMix, is_matrix_like_test) {
  using stan::is_matrix_like;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE((is_matrix_like<Eigen::MatrixXd>::value));
  EXPECT_TRUE((is_matrix_like<Eigen::VectorXd>::value));
  EXPECT_TRUE((is_matrix_like<Eigen::RowVectorXd>::value));
  EXPECT_TRUE((is_matrix_like<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_TRUE((is_matrix_like<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_TRUE((is_matrix_like<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_TRUE((is_matrix_like<var_value<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_matrix_like<var_value<Eigen::VectorXd>>::value));
  EXPECT_TRUE((is_matrix_like<var_value<Eigen::RowVectorXd>>::value));
  Eigen::MatrixXd A(10, 10);
  Eigen::MatrixXd B(10, 10);
  EXPECT_TRUE((is_matrix_like<decltype(A * B)>::value));

  EXPECT_FALSE((is_matrix_like<double>::value));
  EXPECT_FALSE((is_matrix_like<var_value<double>>::value));
}
