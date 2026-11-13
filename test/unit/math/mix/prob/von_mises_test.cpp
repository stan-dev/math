#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

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

TEST_F(AgradRev, mathMixScalFun_von_mises_lpdf_large_kappa) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::von_mises_lpdf;
  // lpdf at y = 2e-6, mu = 0, kappa = 1e12 and its first three
  // kappa-derivatives, from mpmath (60 digits)
  const double y = 2e-6;
  const double kappa = 1e12;
  const double ref[] = {10.896572024760143, -1.4999999999992082e-12,
                        -5.0000000000025e-25, 1.00000000000075e-36};
  auto expect_rel = [](double expected, double actual) {
    EXPECT_NEAR(expected, actual, 1e-13 * stan::math::fabs(expected));
  };

  fvar<fvar<var>> k(fvar<var>(kappa, 1.0), fvar<var>(1.0, 0.0));
  fvar<fvar<var>> lp = von_mises_lpdf(y, 0.0, k);
  expect_rel(ref[0], lp.val_.val_.val());
  expect_rel(ref[1], lp.d_.val_.val());
  expect_rel(ref[1], lp.val_.d_.val());
  expect_rel(ref[2], lp.d_.d_.val());
  lp.d_.d_.grad();
  expect_rel(ref[3], k.val_.val_.adj());

  var kv = kappa;
  var lpv = von_mises_lpdf(y, 0.0, kv);
  lpv.grad();
  expect_rel(ref[1], kv.adj());
}
