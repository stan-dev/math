#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunRev, as_value_array_or_scalar_scalar) {
  double b_val = 4;
  stan::math::var b(b_val);
  EXPECT_EQ(b_val, stan::math::as_value_array_or_scalar(b));
}

TEST(MathFunRev, as_value_array_or_scalar_std_vector_lvalue) {
  int n = 100;
  Eigen::ArrayXd a_val = Eigen::ArrayXd::Random(n);
  std::vector<double> b_val(n);
  for (int i = 0; i < n; i++) {
    b_val[i] = a_val.coeff(i);
  }
  std::vector<stan::math::var> b(b_val.begin(), b_val.end());
  const auto& tmp = stan::math::as_value_array_or_scalar(b);
  Eigen::ArrayXd res = tmp;
  EXPECT_TRUE((stan::is_eigen_array<decltype(tmp)>::value));
  EXPECT_MATRIX_EQ(res, a_val);
}

TEST(MathFunRev, as_value_array_or_scalar_std_vector_rvalue) {
  int n = 100;
  Eigen::ArrayXd a_val = Eigen::ArrayXd::Random(n);
  std::vector<double> b_val(n);
  for (int i = 0; i < n; i++) {
    b_val[i] = a_val.coeff(i);
  }
  std::vector<stan::math::var> b(b_val.begin(), b_val.end());
  const auto& tmp = stan::math::as_value_array_or_scalar(std::move(b));
  EXPECT_TRUE((stan::is_eigen_array<decltype(tmp)>::value));
  Eigen::ArrayXd res = tmp;
  EXPECT_MATRIX_EQ(res, a_val);
}

TEST(MathFunRev, as_value_array_or_scalar_matrix_lvalue) {
  int n = 100;
  Eigen::MatrixXd a_val = Eigen::MatrixXd::Random(n, n);
  Eigen::Matrix<stan::math::var, -1, -1> a(a_val);
  auto&& tmp = stan::math::as_value_array_or_scalar(a);
  EXPECT_TRUE((stan::is_eigen_array<decltype(tmp)>::value));
  Eigen::ArrayXXd res = tmp;
  EXPECT_MATRIX_EQ(res, a_val);
}

TEST(MathFunRev, as_value_array_or_scalar_const_matrix_lvalue) {
  int n = 100;
  const Eigen::MatrixXd a_val = Eigen::MatrixXd::Random(n, n);
  Eigen::Matrix<stan::math::var, -1, -1> a(a_val);
  auto&& tmp = stan::math::as_value_array_or_scalar(a);
  EXPECT_TRUE((stan::is_eigen_array<decltype(tmp)>::value));
  Eigen::ArrayXXd res = tmp;
  EXPECT_MATRIX_EQ(res, a_val);
}

TEST(MathFunRev, as_value_array_or_scalar_matrix_rvalue) {
  int n = 10;
  Eigen::MatrixXd a_val = Eigen::MatrixXd::Random(n, n);
  Eigen::Matrix<stan::math::var, -1, -1> a(a_val);
  auto b = a;
  const auto& tmp = stan::math::as_value_array_or_scalar(std::move(b));
  EXPECT_TRUE((stan::is_eigen_array<decltype(tmp)>::value));
  Eigen::ArrayXXd res = tmp;
  EXPECT_MATRIX_EQ(res, a_val);
}

TEST(MathFunRev, as_value_array_or_scalar_var_value) {
  int n = 100;
  const Eigen::MatrixXd a_val = Eigen::MatrixXd::Random(n, n);
  stan::math::var_value<Eigen::Matrix<double, -1, -1>> a(a_val);
  auto&& tmp = stan::math::as_value_array_or_scalar(a);
  EXPECT_TRUE((stan::is_eigen_array<decltype(tmp)>::value));
  Eigen::ArrayXXd res = tmp;
  EXPECT_MATRIX_EQ(res, a_val);
}
