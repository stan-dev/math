#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, is_var_matrix_test) {
  using stan::is_var_matrix;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_TRUE((is_var_matrix<var_value<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_var_matrix<var_value<Eigen::SparseMatrix<double>>>::value));
  EXPECT_FALSE((is_var_matrix<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_var_matrix<Eigen::Matrix<var_value<float>, -1, -1>>::value));
  EXPECT_FALSE(is_var_matrix<vari>::value);
  EXPECT_FALSE((is_var_matrix<double>::value));
  EXPECT_FALSE((is_var_matrix<vari_value<double>>::value));
  EXPECT_FALSE((is_var_matrix<Eigen::MatrixXd>::value));
}

TEST(MetaTraitsRevScal, is_var_col_vector_test) {
  using stan::is_var_col_vector;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_TRUE(
      (is_var_col_vector<var_value<Eigen::Matrix<double, -1, 1>>>::value));
  EXPECT_FALSE(
      (is_var_col_vector<var_value<Eigen::Matrix<double, 1, -1>>>::value));
  EXPECT_FALSE(
      (is_var_col_vector<var_value<Eigen::Matrix<double, -1, -1>>>::value));
  EXPECT_FALSE((is_var_col_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE(
      (is_var_col_vector<Eigen::Matrix<var_value<float>, -1, 1>>::value));
  EXPECT_FALSE((is_var_col_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_var_col_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_FALSE(
      (is_var_col_vector<Eigen::Matrix<var_value<float>, -1, -1>>::value));
  EXPECT_FALSE(
      (is_var_col_vector<Eigen::Matrix<var_value<float>, 1, -1>>::value));
  EXPECT_FALSE(
      (is_var_col_vector<var_value<Eigen::SparseMatrix<double>>>::value));
  EXPECT_FALSE(is_var_col_vector<vari>::value);
  EXPECT_FALSE((is_var_col_vector<double>::value));
  EXPECT_FALSE((is_var_col_vector<vari_value<double>>::value));
  EXPECT_FALSE((is_var_col_vector<Eigen::MatrixXd>::value));
}

TEST(MetaTraitsRevScal, is_var_row_vector_test) {
  using stan::is_var_row_vector;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_FALSE(
      (is_var_row_vector<var_value<Eigen::Matrix<double, -1, 1>>>::value));
  EXPECT_TRUE(
      (is_var_row_vector<var_value<Eigen::Matrix<double, 1, -1>>>::value));
  EXPECT_FALSE(
      (is_var_row_vector<var_value<Eigen::Matrix<double, -1, -1>>>::value));
  EXPECT_FALSE((is_var_row_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_FALSE(
      (is_var_row_vector<Eigen::Matrix<var_value<float>, 1, -1>>::value));
  EXPECT_FALSE((is_var_row_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE((is_var_row_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE(
      (is_var_row_vector<Eigen::Matrix<var_value<float>, -1, 1>>::value));
  EXPECT_FALSE(
      (is_var_row_vector<Eigen::Matrix<var_value<float>, -1, -1>>::value));
  EXPECT_FALSE(
      (is_var_row_vector<var_value<Eigen::SparseMatrix<double>>>::value));
  EXPECT_FALSE(is_var_row_vector<vari>::value);
  EXPECT_FALSE((is_var_row_vector<double>::value));
  EXPECT_FALSE((is_var_row_vector<vari_value<double>>::value));
  EXPECT_FALSE((is_var_row_vector<Eigen::MatrixXd>::value));
}

TEST(MetaTraitsRevScal, is_var_vector_test) {
  using stan::is_var_vector;
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::vari;
  using stan::math::vari_value;
  EXPECT_TRUE((is_var_vector<var_value<Eigen::Matrix<double, -1, 1>>>::value));
  EXPECT_TRUE((is_var_vector<var_value<Eigen::Matrix<double, 1, -1>>>::value));
  EXPECT_FALSE(
      (is_var_vector<var_value<Eigen::Matrix<double, -1, -1>>>::value));
  EXPECT_FALSE((is_var_vector<Eigen::Matrix<var, 1, -1>>::value));
  EXPECT_FALSE((is_var_vector<Eigen::Matrix<var, -1, 1>>::value));
  EXPECT_FALSE((is_var_vector<Eigen::Matrix<var, -1, -1>>::value));
  EXPECT_FALSE((is_var_vector<var_value<Eigen::SparseMatrix<double>>>::value));
  EXPECT_FALSE(is_var_vector<vari>::value);
  EXPECT_FALSE((is_var_vector<double>::value));
  EXPECT_FALSE((is_var_vector<vari_value<double>>::value));
  EXPECT_FALSE((is_var_vector<Eigen::MatrixXd>::value));
}
