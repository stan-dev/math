#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <gtest/gtest.h>

TEST(MathMetaMix, is_var_eigen_test) {
  using stan::is_var_eigen;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE((is_var_eigen<var_value<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_var_eigen<var_value<Eigen::ArrayXd>>::value));
  EXPECT_TRUE((is_var_eigen<var_value<Eigen::VectorXd>>::value));
  EXPECT_TRUE((is_var_eigen<var_value<Eigen::RowVectorXd>>::value));

  Eigen::MatrixXd A(10, 10);
  Eigen::MatrixXd B(10, 10);
  EXPECT_FALSE((is_var_eigen<decltype(A * B)>::value));
  EXPECT_FALSE((is_var_eigen<Eigen::MatrixXd>::value));
  EXPECT_FALSE((is_var_eigen<Eigen::VectorXd>::value));
  EXPECT_FALSE((is_var_eigen<Eigen::RowVectorXd>::value));
  EXPECT_FALSE((is_var_eigen<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_var_eigen<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE((is_var_eigen<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_FALSE((is_var_eigen<double>::value));
  EXPECT_FALSE((is_var_eigen<var_value<double>>::value));
}
