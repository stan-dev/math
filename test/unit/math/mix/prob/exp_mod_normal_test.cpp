#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_lpdf) {
  auto f = [](const auto& p) {
    return stan::math::exp_mod_normal_lpdf(p[0], p[1], p[2], p[3]);
  };

  for (double y : {-40.0, -10.0, -3.0, 0.3, 4.0, 10.0, 20.0, 39.0, 60.0}) {
    stan::test::expect_ad(f, Eigen::Vector4d(y, 0.1, 1.3, 0.7));
  }
  stan::test::expect_ad(f, Eigen::Vector4d(0, 0, 2, 20));
  stan::test::expect_ad(
      [](const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lpdf(Eigen::Vector3d(-40, 0.3, 60),
                                               mu, sigma, lambda);
      },
      0.1, 1.3, 0.7);
}

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_lpdf_vector_tails) {
  // z = (y - mu) / sigma - lambda sigma runs from the deep lower tail to
  // the positive tail.
  Eigen::VectorXd lambda(5);
  lambda << 0.05, 0.7, 2.5, 20, 40;
  stan::test::expect_ad(
      [](const auto& y, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lpdf(y, 0.1, sigma, lambda);
      },
      0.3, 1.3, lambda);
  Eigen::VectorXd mu(4);
  mu << -40, 0, 10, 39;
  Eigen::VectorXd sigma(4);
  sigma << 0.5, 1, 2, 3;
  stan::test::expect_ad(
      [](const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lpdf(0.3, mu, sigma, lambda);
      },
      mu, sigma, 0.7);
}

TEST_F(AgradRev, mathMixScalFun_exp_mod_normal_lpdf_large_lambda) {
  // MPFR reference; erfc underflows here.
  EXPECT_NEAR(-1.6127097401997972,
              stan::math::exp_mod_normal_lpdf(0.0, 0.0, 2.0, 20.0), 1e-12);
}
