#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, as_array_or_scalar_scalar) {
  int a = 3;
  double b = 4;
  EXPECT_EQ(a, stan::math::as_array_or_scalar(a));
  EXPECT_EQ(b, stan::math::as_array_or_scalar(b));
}

TEST(MathMetaPrim, as_array_or_scalar_std_vector_lvalue) {
  int n = 100;
  Eigen::ArrayXd a = Eigen::ArrayXd::Random(n);
  std::vector<double> b(n);
  for (int i = 0; i < n; i++) {
    b[i] = a.coeff(i);
  }
  const auto& tmp = stan::math::as_array_or_scalar(b);
  Eigen::ArrayXd res = tmp;
  EXPECT_MATRIX_EQ(res, a);
}

TEST(MathMetaPrim, as_array_or_scalar_std_vector_rvalue) {
  int n = 100;
  Eigen::ArrayXd a = Eigen::ArrayXd::Random(n);
  std::vector<double> b(n);
  for (int i = 0; i < n; i++) {
    b[i] = a.coeff(i);
  }
  const auto& tmp = stan::math::as_array_or_scalar(std::move(b));
  b.assign(n, 0);  // overwrite the memory if b was not moved
  Eigen::ArrayXd res = tmp;
  EXPECT_MATRIX_EQ(res, a);
}

TEST(MathMetaPrim, as_array_or_scalar_matrix_lvalue) {
  int n = 100;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  auto&& tmp = stan::math::as_array_or_scalar(a);
  tmp(0, 0) = 1234;
  Eigen::ArrayXXd res = tmp;
  EXPECT_MATRIX_EQ(res, a);
}

TEST(MathMetaPrim, as_array_or_scalar_const_matrix_lvalue) {
  int n = 100;
  const Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  const auto& tmp = stan::math::as_array_or_scalar(a);
  Eigen::ArrayXXd res = tmp;
  EXPECT_MATRIX_EQ(res, a);
}

TEST(MathMetaPrim, as_array_or_scalar_matrix_rvalue) {
  int n = 10;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd b = a;
  const auto& tmp = stan::math::as_array_or_scalar(std::move(b));
  b.setZero();  // overwrite the memory if b was not moved
  Eigen::ArrayXXd res = tmp;
  EXPECT_MATRIX_EQ(res, a);
}
