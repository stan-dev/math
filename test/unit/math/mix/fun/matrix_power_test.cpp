#include <stan/math/rev.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrixPower, ad_tests) {
  using Eigen::MatrixXd;
  using stan::math::matrix_power;
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto fn2 = [](const auto& M) { return matrix_power(M, -2); };
  auto f0 = [](const auto& M) { return matrix_power(M, 0); };
  auto f2 = [](const auto& M) { return matrix_power(M, 2); };
  auto f1 = [](const auto& M) { return matrix_power(M, 1); };
  auto f3 = [](const auto& M) { return matrix_power(M, 3); };
  auto f9 = [](const auto& M) { return matrix_power(M, 9); };

  MatrixXd M11(1, 1);
  M11 << -5.0;
  expect_ad(f0, M11);
  expect_ad(f1, M11);
  expect_ad(f2, M11);
  expect_ad(f3, M11);
  expect_ad(f9, M11);

  expect_ad_matvar(f0, M11);
  expect_ad_matvar(f1, M11);
  expect_ad_matvar(f2, M11);
  expect_ad_matvar(f3, M11);
  expect_ad_matvar(f9, M11);

  MatrixXd M22(2, 2);
  M22 << 10.0, -2.0, 3.5, 4.0;
  expect_ad(f0, M22);
  expect_ad(f1, M22);
  expect_ad(f2, M22);
  expect_ad(f3, M22);
  expect_ad(f9, M22);

  expect_ad_matvar(f0, M22);
  expect_ad_matvar(f1, M22);
  expect_ad_matvar(f2, M22);
  expect_ad_matvar(f3, M22);
  expect_ad_matvar(f9, M22);

  MatrixXd M33(3, 3);
  M33 << 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, -0.5, 2.0, 3.0;
  expect_ad(f0, M33);
  expect_ad(f1, M33);
  expect_ad(f2, M33);
  expect_ad(f3, M33);

  expect_ad_matvar(f0, M33);
  expect_ad_matvar(f1, M33);
  expect_ad_matvar(f2, M33);
  expect_ad_matvar(f3, M33);
  // Finite differences too inaccurate.
  // expect_ad(f9, M33);

  expect_ad(fn2, M33);
  expect_ad_matvar(fn2, M33);

  MatrixXd not_square = Eigen::MatrixXd::Identity(3, 4);
  expect_ad(f2, not_square);
  expect_ad_matvar(f2, not_square);

  Eigen::MatrixXd zero_size = Eigen::MatrixXd::Identity(0, 0);
  expect_ad(f2, zero_size);
  expect_ad_matvar(f2, zero_size);

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -inf;

  MatrixXd nan_not_square = not_square;
  nan_not_square(0, 0) = nan;
  expect_ad(fn2, nan_not_square);
  expect_ad_matvar(fn2, nan_not_square);

  MatrixXd nan_33 = M33;
  nan_33(0, 0) = nan;
  expect_ad(fn2, nan_33);
  expect_ad(f2, nan_33);
  expect_ad_matvar(fn2, nan_33);
  expect_ad_matvar(f2, nan_33);

  MatrixXd inf_33 = M33;
  inf_33(0, 0) = inf;
  expect_ad(f2, inf_33);
  expect_ad_matvar(f2, inf_33);

  MatrixXd ninf_33 = M33;
  ninf_33(0, 0) = ninf;
  expect_ad(f2, ninf_33);
  expect_ad_matvar(f2, ninf_33);

  stan::math::recover_memory();
}
