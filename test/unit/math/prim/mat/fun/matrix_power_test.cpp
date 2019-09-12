#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPower, matrix_power_two_by_two) {
  using stan::math::matrix_power;

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
  Eigen::MatrixXd M(2, 2);
  M << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd M2(2, 2);
  M2 << 7.0, 10.0, 15.0, 22.0;
  Eigen::MatrixXd M3(2, 2);
  M3 << 37.0, 54.0, 81.0, 118.0;

  expect_matrix_eq(I, matrix_power(M, 0));
  expect_matrix_eq(M, matrix_power(M, 1));
  expect_matrix_eq(M2, matrix_power(M, 2));
  expect_matrix_eq(M3, matrix_power(M, 3));
}

TEST(MathMatrixPower, matrix_power_one_by_one) {
  using stan::math::matrix_power;

  Eigen::MatrixXd I(1, 1);
  I << 1.0;
  Eigen::MatrixXd M(1, 1);
  M << 3.0;
  Eigen::MatrixXd M2(1, 1);
  M2 << 9.0;
  Eigen::MatrixXd M3(1, 1);
  M3 << 27.0;

  expect_matrix_eq(I, matrix_power(M, 0));
  expect_matrix_eq(M, matrix_power(M, 1));
  expect_matrix_eq(M2, matrix_power(M, 2));
  expect_matrix_eq(M3, matrix_power(M, 3));
}

TEST(MathMatrixPower, matrix_power_compare_to_simple_impl) {
  using stan::math::matrix_power;

  int size = 3;
  Eigen::MatrixXd M = Eigen::MatrixXd::Random(size, size);
  Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(size, size);
  int exponent = 3;
  for (int i = 0; i < exponent; i++)
    expected *= M;

  expect_matrix_eq(expected, matrix_power(M, exponent));
}

TEST(MathMatrixPower, matrix_power_invalid) {
  using stan::math::matrix_power;

  Eigen::MatrixXd not_square = Eigen::MatrixXd::Identity(3, 4);
  Eigen::MatrixXd good = Eigen::MatrixXd::Identity(2, 2);

  EXPECT_THROW(matrix_power(not_square, 2), std::invalid_argument);
  EXPECT_THROW(matrix_power(good, -2), std::invalid_argument);
}
