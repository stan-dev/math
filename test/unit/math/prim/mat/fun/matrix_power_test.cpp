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

TEST(MathMatrixPower, matrix_power_prim_one_by_one) {
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

TEST(MathMatrixPower, matrix_power_prim_compare_to_simple_impl) {
  using stan::math::matrix_power;

  int size = 5;
  Eigen::MatrixXd M = Eigen::MatrixXd::Random(size, size);
  Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(size, size);
  int exponent = 4;
  for (int i = 0; i < exponent; i++)
    expected *= M;

  expect_matrix_eq(expected, matrix_power(M, exponent));
}

TEST(MathMatrixPower, matrix_power_prim_invalid) {
  using stan::math::matrix_power;

  Eigen::MatrixXd zero_size = Eigen::MatrixXd::Identity(0, 0);
  Eigen::MatrixXd not_square = Eigen::MatrixXd::Identity(3, 4);
  Eigen::MatrixXd good = Eigen::MatrixXd::Identity(2, 2);

  EXPECT_NO_THROW(matrix_power(good, 2));

  EXPECT_THROW(matrix_power(zero_size, 2), std::invalid_argument);
  EXPECT_THROW(matrix_power(not_square, 2), std::invalid_argument);
  EXPECT_THROW(matrix_power(good, -2), std::invalid_argument);

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -inf;

  // invalid_argument takes precedence over domain_error.
  not_square(0, 0) = nan;
  good(0, 0) = nan;
  EXPECT_THROW(matrix_power(not_square, 2), std::invalid_argument);
  EXPECT_THROW(matrix_power(good, -2), std::invalid_argument);

  good(0, 0) = nan;
  EXPECT_THROW(matrix_power(good, 2), std::domain_error);
  good(0, 0) = inf;
  EXPECT_THROW(matrix_power(good, 2), std::domain_error);
  good(0, 0) = ninf;
  EXPECT_THROW(matrix_power(good, 2), std::domain_error);
}
