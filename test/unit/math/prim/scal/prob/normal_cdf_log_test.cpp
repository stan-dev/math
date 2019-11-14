#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(ProbNormal, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::normal_lcdf(y, mu, sigma)),
                  (stan::math::normal_cdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::normal_lcdf<double, double, double>(y, mu, sigma)),
      (stan::math::normal_cdf_log<double, double, double>(y, mu, sigma)));
}

TEST(mathMixScalFun, lcdf_derivatives) {
  auto f = [](const double mu, const double sigma) {
    return [=](const auto& y) { return stan::math::normal_lcdf(y, mu, sigma); };
  };
  stan::test::expect_ad(f(0.0, 1.0), -50.0);
  stan::test::expect_ad(f(0.0, 1.0), -20.0 * stan::math::SQRT_2);
  stan::test::expect_ad(f(0.0, 1.0), -5.5);
  stan::test::expect_ad(f(0.0, 1.0), 0.0);
  stan::test::expect_ad(f(0.0, 1.0), 0.15);
  stan::test::expect_ad(f(0.0, 1.0), 1.14);
  stan::test::expect_ad(f(0.0, 1.0), 3.00);
  stan::test::expect_ad(f(0.0, 1.0), 10.00);
  stan::test::expect_ad(f(-1.0, 2.0), 1.50);
  stan::test::expect_ad(f(2.0, 1.0), 0.50);
}
