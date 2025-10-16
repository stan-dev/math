#include <stan/math/rev.hpp>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctionsPromoteScalar, VarMatrix) {
  using stan::math::var;
  using stan::math::var_value;
  expect_promote_type<var_value<Eigen::MatrixXd>, var,
                      var_value<Eigen::MatrixXd>>();
  expect_promote_type<var_value<Eigen::VectorXd>, var,
                      var_value<Eigen::VectorXd>>();
  expect_promote_type<var_value<Eigen::RowVectorXd>, var,
                      var_value<Eigen::RowVectorXd>>();
  expect_promote_type<Eigen::MatrixXd, double, var_value<Eigen::MatrixXd>>();
  expect_promote_type<Eigen::VectorXd, double, var_value<Eigen::VectorXd>>();
  expect_promote_type<Eigen::RowVectorXd, double,
                      var_value<Eigen::RowVectorXd>>();
}
