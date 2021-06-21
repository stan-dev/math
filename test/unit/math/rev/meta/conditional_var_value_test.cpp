#include <stan/math/rev/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaRev, conditional_var_value_scalar) {
  using stan::conditional_var_value_t;
  using stan::math::var;
  EXPECT_SAME_TYPE(double, conditional_var_value_t<double, double>);
  EXPECT_SAME_TYPE(var, conditional_var_value_t<var, double>);
}

TEST(MathMetaRev, conditional_var_value_reference) {
  using stan::conditional_var_value_t;
  EXPECT_SAME_TYPE(double, conditional_var_value_t<double, double>);
  EXPECT_SAME_TYPE(double, conditional_var_value_t<double&, double&>);
  EXPECT_SAME_TYPE(double, conditional_var_value_t<double&&, double&&>);
}

TEST(MathMetaRev, conditional_var_value_vector) {
  using stan::conditional_var_value_t;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_SAME_TYPE(Eigen::VectorXd,
                   conditional_var_value_t<double, Eigen::VectorXd>);
  EXPECT_SAME_TYPE(var_value<Eigen::VectorXd>,
                   conditional_var_value_t<var, Eigen::VectorXd>);
}

TEST(MathMetaRev, conditional_var_value_row_vector) {
  using stan::conditional_var_value_t;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_SAME_TYPE(Eigen::RowVectorXd,
                   conditional_var_value_t<double, Eigen::RowVectorXd>);
  EXPECT_SAME_TYPE(var_value<Eigen::RowVectorXd>,
                   conditional_var_value_t<var, Eigen::RowVectorXd>);
}

TEST(MathMetaRev, conditional_var_value_matrix) {
  using stan::conditional_var_value_t;
  using stan::math::var;
  using stan::math::var_value;
  EXPECT_SAME_TYPE(Eigen::MatrixXd,
                   conditional_var_value_t<double, Eigen::MatrixXd>);
  EXPECT_SAME_TYPE(var_value<Eigen::MatrixXd>,
                   conditional_var_value_t<var, Eigen::MatrixXd>);
  EXPECT_SAME_TYPE(var_value<Eigen::MatrixXd>,
                   conditional_var_value_t<var, Eigen::Matrix<var, -1, -1>>);
}

TEST(MathMetaRev, conditional_var_value_expression) {
  using stan::conditional_var_value_t;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd a;
  Eigen::MatrixXd b;
  EXPECT_SAME_TYPE(
      Eigen::MatrixXd,
      conditional_var_value_t<double, decltype(a + b.block(1, 1, 2, 2))>);
  EXPECT_SAME_TYPE(
      var_value<Eigen::MatrixXd>,
      conditional_var_value_t<var, decltype(a + b.block(1, 1, 2, 2))>);
}

TEST(MathMetaRev, conditional_var_value_container_T_scalar) {
  using stan::conditional_var_value_t;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd a;
  Eigen::MatrixXd b;
  EXPECT_SAME_TYPE(Eigen::MatrixXd,
                   conditional_var_value_t<Eigen::VectorXd, Eigen::MatrixXd>);
  EXPECT_SAME_TYPE(
      var_value<Eigen::MatrixXd>,
      conditional_var_value_t<Eigen::Matrix<var, 1, Eigen::Dynamic>,
                              Eigen::MatrixXd>);
}
