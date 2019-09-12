#include <stan/math/rev/mat.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPower, matrix_power_ad_tests) {
  using Eigen::MatrixXd;
  using stan::math::matrix_power;
  using stan::test::expect_ad;

  auto f0 = [](const auto& M) { return matrix_power(M, 0); };
  auto f1 = [](const auto& M) { return matrix_power(M, 1); };
  auto f3 = [](const auto& M) { return matrix_power(M, 3); };

  MatrixXd M11(1, 1);
  M11 << -5.0;
  expect_ad(f0, M11);
  expect_ad(f1, M11);
  expect_ad(f3, M11);

  MatrixXd M22(2, 2);
  M22 << 10.0, -2.0, 3.5, 4.0;
  expect_ad(f0, M22);
  expect_ad(f1, M22);
  expect_ad(f3, M22);

  MatrixXd M33(3, 3);
  M33 << 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, -0.5, 2.0, 3.0;
  expect_ad(f0, M33);
  expect_ad(f1, M33);
  expect_ad(f3, M33);

  MatrixXd not_square = Eigen::MatrixXd::Identity(3, 4);
  expect_ad(f0, not_square);
  expect_ad(f1, not_square);
  expect_ad(f3, not_square);

  auto fn2 = [](const auto& M) { return matrix_power(M, -2); };
  expect_ad(fn2, M11);

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -inf;

  Eigen::MatrixXd zero_size = Eigen::MatrixXd::Identity(0, 0);
  expect_ad(f3, zero_size);

  // invalid_argument takes precedence over domain_error.
  not_square(0, 0) = nan;
  expect_ad(f3, not_square);
  M33(0, 0) = nan;
  expect_ad(fn2, M33);

  M33(0, 0) = nan;
  expect_ad(f3, M33);
  M33(0, 0) = inf;
  expect_ad(f3, M33);
  M33(0, 0) = ninf;
  expect_ad(f3, M33);
}
