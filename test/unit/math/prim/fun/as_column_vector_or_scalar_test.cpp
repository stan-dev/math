#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, as_column_vector_or_scalar_scalar) {
  using stan::math::as_column_vector_or_scalar;
  int a = 3;
  double b = 4;
  EXPECT_EQ(a, as_column_vector_or_scalar(a));
  EXPECT_EQ(b, as_column_vector_or_scalar(b));
}

TEST(MathMetaPrim, as_column_vector_or_scalar_std_vector_lvalue) {
  using stan::math::as_column_vector_or_scalar;
  int n = 100;
  Eigen::VectorXd a = Eigen::VectorXd::Random(n);
  std::vector<double> b(n);
  for (int i = 0; i < n; i++) {
    b[i] = a.coeff(i);
  }
  auto&& tmp = as_column_vector_or_scalar(b);
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp)>::value));
  a[0] = tmp[0] = 12345;
  Eigen::VectorXd res = tmp;
  EXPECT_MATRIX_EQ(res, a);
}

TEST(MathMetaPrim, as_column_vector_or_scalar_std_vector_rvalue) {
  using stan::math::as_column_vector_or_scalar;
  int n = 100;
  Eigen::VectorXd a = Eigen::VectorXd::Random(n);
  std::vector<double> b(n);
  for (int i = 0; i < n; i++) {
    b[i] = a.coeff(i);
  }
  const auto& tmp = as_column_vector_or_scalar(std::move(b));
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp)>::value));
  b.assign(n, 0);  // overwrite the memory if b was not moved
  Eigen::VectorXd res = tmp;
  EXPECT_MATRIX_EQ(res, a);
}

TEST(MathMetaPrim, as_column_vector_or_scalar_const_row_vector_lvalue) {
  using stan::math::as_column_vector_or_scalar;
  int n = 100;
  const Eigen::RowVectorXd a = Eigen::RowVectorXd::Random(n);
  auto&& tmp = as_column_vector_or_scalar(a);
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp)>::value));
  Eigen::VectorXd res = tmp;
  Eigen::VectorXd correct = a.transpose();
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(MathMetaPrim, as_column_vector_or_scalar_row_vector_lvalue) {
  using stan::math::as_column_vector_or_scalar;
  int n = 100;
  Eigen::RowVectorXd a = Eigen::RowVectorXd::Random(n);
  auto&& tmp = as_column_vector_or_scalar(a);
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp)>::value));
  tmp[0] = 1234;
  Eigen::VectorXd res = tmp;
  Eigen::VectorXd correct = a.transpose();
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(MathMetaPrim, as_column_vector_or_scalar_row_vector_rvalue) {
  using stan::math::as_column_vector_or_scalar;
  int n = 100;
  Eigen::RowVectorXd a = Eigen::RowVectorXd::Random(n);
  Eigen::RowVectorXd b = a;
  const auto& tmp = as_column_vector_or_scalar(std::move(b));
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp)>::value));
  Eigen::VectorXd res = tmp;
  b.setZero();  // overwrite the memory if b was not moved
  Eigen::VectorXd correct = a.transpose();
  EXPECT_MATRIX_EQ(res, correct);
}
