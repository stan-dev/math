#include <stan/math/rev/mat.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPower, matrix_power_ad_tests) {
  using stan::math::matrix_power;
  auto f0 = [](const auto& M) { return matrix_power(M, 0); };
  auto f1 = [](const auto& M) { return matrix_power(M, 1); };
  auto f3 = [](const auto& M) { return matrix_power(M, 3); };
  Eigen::MatrixXd M11(1, 1);
  M11 << -5.0;
  stan::test::expect_ad(f0, M11);
  stan::test::expect_ad(f1, M11);
  stan::test::expect_ad(f3, M11);
  Eigen::MatrixXd M22(2, 2);
  M22 << 10.0, -2.0, 3.5, 4.0;
  stan::test::expect_ad(f0, M22);
  stan::test::expect_ad(f1, M22);
  stan::test::expect_ad(f3, M22);
  Eigen::MatrixXd M33(3, 3);
  M33 << 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, -0.5, 2.0, 3.0;
  stan::test::expect_ad(f0, M33);
  stan::test::expect_ad(f1, M33);
  stan::test::expect_ad(f3, M33);
}
