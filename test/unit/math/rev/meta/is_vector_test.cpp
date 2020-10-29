#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, is_col_vector_test) {
  using stan::is_col_vector;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_TRUE((is_col_vector<var_value<Eigen::Matrix<double, -1, 1>>>::value));
  EXPECT_TRUE((is_col_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_TRUE((is_col_vector<Eigen::Matrix<var_value<float>, -1, 1>>::value));
  EXPECT_FALSE((is_col_vector<var_value<Eigen::Matrix<double, 1, -1>>>::value));
  EXPECT_FALSE(
      (is_col_vector<var_value<Eigen::Matrix<double, -1, -1>>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::Matrix<var_value<float>, -1, -1>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::Matrix<var_value<float>, 1, -1>>::value));
  EXPECT_FALSE((is_col_vector<var_value<Eigen::SparseMatrix<double>>>::value));
  EXPECT_FALSE(is_col_vector<vari>::value);
  EXPECT_FALSE((is_col_vector<double>::value));
  EXPECT_FALSE((is_col_vector<vari_value<double>>::value));
  EXPECT_FALSE((is_col_vector<Eigen::MatrixXd>::value));
}

TEST(MetaTraitsRevScal, is_row_vector_test) {
  using stan::is_row_vector;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_TRUE((is_row_vector<var_value<Eigen::Matrix<double, 1, -1>>>::value));
  EXPECT_TRUE((is_row_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_TRUE((is_row_vector<Eigen::Matrix<var_value<float>, 1, -1>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::Matrix<var_value<float>, -1, 1>>::value));
  EXPECT_FALSE((is_row_vector<var_value<Eigen::Matrix<double, -1, 1>>>::value));
  EXPECT_FALSE(
      (is_row_vector<var_value<Eigen::Matrix<double, -1, -1>>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::Matrix<var_value<float>, -1, -1>>::value));
  EXPECT_FALSE((is_row_vector<var_value<Eigen::SparseMatrix<double>>>::value));
  EXPECT_FALSE(is_row_vector<vari>::value);
  EXPECT_FALSE((is_row_vector<double>::value));
  EXPECT_FALSE((is_row_vector<vari_value<double>>::value));
  EXPECT_FALSE((is_row_vector<Eigen::MatrixXd>::value));
}

TEST(MetaTraitsRevScal, is_vector_test) {
  using stan::is_vector;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_TRUE((is_vector<var_value<Eigen::Matrix<double, -1, 1>>>::value));
  EXPECT_TRUE((is_vector<var_value<Eigen::Matrix<double, 1, -1>>>::value));
  EXPECT_TRUE((is_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_TRUE((is_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE((is_vector<var_value<Eigen::Matrix<double, -1, -1>>>::value));
  EXPECT_FALSE((is_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_vector<var_value<Eigen::SparseMatrix<double>>>::value));
  EXPECT_FALSE(is_vector<vari>::value);
  EXPECT_FALSE((is_vector<double>::value));
  EXPECT_FALSE((is_vector<vari_value<double>>::value));
  EXPECT_FALSE((is_vector<Eigen::MatrixXd>::value));
}
