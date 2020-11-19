#include <stan/math/mix.hpp>
#include <test/unit/math/mix/meta/eigen_utils.hpp>
#include <gtest/gtest.h>

TEST(MathMetaMix, is_var_dense_dynamic_test) {
  using stan::is_var_dense_dynamic;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_TRUE((is_var_dense_dynamic<var_value<Eigen::MatrixXd>>::value));
  EXPECT_TRUE(
      (is_var_dense_dynamic<var_value<Eigen::Array<double, -1, -1>>>::value));
  EXPECT_FALSE((is_var_dense_dynamic<var_value<Eigen::VectorXd>>::value));
  EXPECT_FALSE((is_var_dense_dynamic<var_value<Eigen::RowVectorXd>>::value));

  Eigen::MatrixXd A(10, 10);
  Eigen::MatrixXd B(10, 10);
  EXPECT_FALSE((is_var_dense_dynamic<decltype(A * B)>::value));
  EXPECT_FALSE((is_var_dense_dynamic<Eigen::MatrixXd>::value));
  EXPECT_FALSE((is_var_dense_dynamic<Eigen::VectorXd>::value));
  EXPECT_FALSE((is_var_dense_dynamic<Eigen::RowVectorXd>::value));
  EXPECT_FALSE((is_var_dense_dynamic<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_var_dense_dynamic<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE((is_var_dense_dynamic<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_FALSE((is_var_dense_dynamic<double>::value));
  EXPECT_FALSE((is_var_dense_dynamic<var_value<double>>::value));
}
