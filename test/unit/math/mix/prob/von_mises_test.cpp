#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST_F(AgradRev, mathMixScalFun_von_mises_lpdf) {
  auto f = [](const auto& y, const auto& mu, const auto& kappa) {
    return stan::math::von_mises_lpdf(y, mu, kappa);
  };

  Eigen::VectorXd y(2);
  y << 1.0, 2.0;
  Eigen::VectorXd mu(2);
  mu << 1.0, 0.5;
  Eigen::VectorXd kappa(2);
  kappa << 1.0, 0.5;

  stan::test::expect_ad(f, y[0], mu[0], kappa[0]);
  stan::test::expect_ad(f, y[0], mu, kappa);
  stan::test::expect_ad(f, y, mu[0], kappa);
  stan::test::expect_ad(f, y, mu, kappa[0]);
  stan::test::expect_ad(f, y[0], mu[0], kappa);
  stan::test::expect_ad(f, y, mu[0], kappa[0]);
  stan::test::expect_ad(f, y[0], mu, kappa[0]);
}
